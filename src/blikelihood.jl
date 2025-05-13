
function b_unpack_params(theta, dimvals, ranks)
    p_dimvals = prod(dimvals)
    num_delta = ranks[1] * (dimvals[1] - ranks[1])
    num_gamma = ranks[2] * (dimvals[2] - ranks[2])
    num_u3 = ranks[1] * dimvals[1] - 1
    num_u4 = ranks[2] * dimvals[2]
    num_sigma = Int(p_dimvals * (p_dimvals + 1) / 2)

    @assert length(theta) == (num_delta + num_gamma + num_u3 + num_u4 + num_sigma) "Parameter vector has wrong length!"

    idx_gamma = num_delta + num_gamma
    idx_u3 = num_delta + num_gamma + num_u3
    idx_u4 = num_delta + num_gamma + num_u3 + num_u4

    delta_vec = theta[1:num_delta]
    gamma_vec = theta[(num_delta+1):idx_gamma]
    u3_vec = theta[(idx_gamma+1):idx_u3]
    u4_vec = theta[(idx_u3+1):idx_u4]
    ll_vec = theta[(idx_u4+1):end]

    delta_star = reshape(delta_vec, ranks[1], dimvals[1] - ranks[1])
    delta = vcat(I(dimvals[1] - ranks[1]), delta_star)

    gamma_star = reshape(gamma_vec, ranks[2], dimvals[2] - ranks[2])
    gamma = vcat(I(dimvals[2] - ranks[2]), gamma_star)

    insert!(u3_vec, 1, 1)
    u3 = reshape(u3_vec, dimvals[1], ranks[1])
    u4 = reshape(u4_vec, dimvals[2], ranks[2])

    ll = vec_to_ll(ll_vec, p_dimvals)

    return (; delta, gamma, u3, u4, ll)
end

function b_pack_params(delta, gamma, u3, u4, ll)
    u3_removed = vec(u3)[2:end]
    return vcat(vec(delta), vec(gamma), u3_removed, vec(u4), vech(ll))
end

function rand_init(dimvals, ranks; spread=1)
    num_delta = ranks[1] * (dimvals[1] - ranks[1])
    num_gamma = ranks[2] * (dimvals[2] - ranks[2])
    num_u3 = ranks[1] * dimvals[1]
    num_u4 = ranks[2] * dimvals[2]

    delta_init = spread .* randn(num_delta)
    gamma_init = spread .* randn(num_gamma)
    u3_init = spread .* randn(num_u3)
    s = u3_init[1]
    u3_new = u3_init / s
    u4_init = spread .* randn(num_u4) * s

    return b_pack_params(delta_init, gamma_init, u3_new, u4_init, I(prod(dimvals)))

end

function init_both(resp, pred, dimvals, ranks; pack_params=true)

    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    coef = ols_coef(resp, pred)

    tensor_coef =
        matten(coef, [1, 2], [3, 4], [dimvals[1], dimvals[2], dimvals[1], dimvals[2]])
    flat1 = tenmat(tensor_coef, row=1)
    flat2 = tenmat(tensor_coef, row=2)
    flat3 = tenmat(tensor_coef, row=3)
    flat4 = tenmat(tensor_coef, row=4)
    u1 = svd(flat1).U[:, 1:ranks[1]]
    u2 = svd(flat2).U[:, 1:ranks[2]]
    u3 = svd(flat3).U[:, 1:ranks[1]]
    u4 = svd(flat4).U[:, 1:ranks[2]]

    s = u3[1, 1]
    u3 = copy(u3) / s
    u4 = copy(u4) * s

    delta = nullspace(u1')
    delta_rot = delta * inv(delta[1:(N1_r1), 1:N1_r1])
    delta_star = delta_rot[(N1_r1+1):end, :]

    gamma = nullspace(u2')
    gamma_rot = gamma * inv(gamma[1:(N2_r2), 1:N2_r2])
    gamma_star = gamma_rot[(N2_r2+1):end, :]

    if pack_params
        return b_pack_params(delta_star, gamma_star, u3, u4, I(prod(dimvals)))
    else
        return (; u1, u2, u3, u4)
    end

end

function both_loglike(theta, resp, pred, dimvals, ranks)
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    delta_rot, gamma_rot, u3, u4, ll = b_unpack_params(theta, dimvals, ranks)
    delta_star = delta_rot[(N1_r1+1):end, :]
    gamma_star = gamma_rot[(N2_r2+1):end, :]

    pi_mat = pi_from_both(u3, u4, dimvals, ranks)

    obs = size(resp, 2)
    omega = omega_from_both(delta_star, gamma_star, dimvals, ranks)
    sparse_omega = sparse(omega)
    sparse_pi = sparse(pi_mat)
    det_term = det(0.5 .* (ll * ll' + ll * ll'))
    if det_term <= 0.0
        return 1e10
    end

    logdet_term = log(det_term)
    X = sparse_omega * ll
    precision_matrix = inv(X') * inv(X)

    sse = 0.0

    for i = 2:obs
        yt = resp[:, i]
        yt_m1 = pred[:, i]
        resid = sparse_omega * yt - sparse_pi * yt_m1
        sse += dot(resid, precision_matrix * resid)
    end

    return 0.5 * ((obs - 1) * logdet_term + sse)
end

function both_hess(theta_est, resp, pred, dimvals, ranks)
    grad_est = ForwardDiff.hessian(
        theta -> both_loglike(theta, resp, pred, dimvals, ranks),
        theta_est,
    )
    return grad_est
end

function comovement_init(resp, pred, dimvals, ranks; iters=3, tol=1e-01, num_starts=20, num_selected=5, spread=1)
    some_init = init_both(resp, pred, dimvals, ranks)
    init_length = length(some_init)
    potential_starts = fill(NaN, init_length + 1, num_starts)
    obj = tet -> both_loglike(tet, resp, pred, dimvals, ranks)
    num_iters = zeros(num_starts)
    problem_starts = zeros(num_starts)

    for i in 1:num_starts
        if i == 1
            both_init = copy(some_init)
        else
            #=both_init = copy(some_init) .+ rand_init(dimvals, ranks)=#
            both_init = rand_init(dimvals, ranks; spread)
        end
        potential_starts[1:(end-1), i] = both_init
        td = TwiceDifferentiable(obj, both_init, autodiff=:forward)

        res = optimize(
            td,
            both_init,
            #=GradientDescent(),=#
            NewtonTrustRegion(;
                initial_delta=1e2,   # Δ₀
                delta_hat=1e4,  # max Δₖ
                eta=0.01,    # accept step when ρₖ > η
                rho_lower=0.1,   # shrink when ρₖ < ρ_lower
                rho_upper=0.9,   # grow   when ρₖ > ρ_upper
            ),
            Optim.Options(iterations=iters, f_abstol=tol, f_reltol=tol, g_abstol=1.0),
        )
        if res.g_residual > 1.0
            problem_starts[i] = 1
        end
        num_iters[i] = res.iterations
        potential_starts[end, i] = res.minimum
        potential_starts[1:(end-1), i] = res.minimizer
    end
    chosen_idx = partialsortperm(potential_starts[end, :], 1:num_selected)
    chosen_start = potential_starts[1:(end-1), chosen_idx]

    return (; chosen_start, num_iters, problem_starts)

end

function main_algorithm(resp, pred, dimvals, ranks; iters=100, tol=1e-05, num_starts=20, num_selected=5, spread=1)

    obj = tet -> both_loglike(tet, resp, pred, dimvals, ranks)
    td = nothing
    res = nothing

    chosen_start, num_iters, problem_starts = comovement_init(resp, pred, dimvals, ranks; iters=3, tol=1e-01, num_starts=num_starts, num_selected=num_selected, spread)
    potential_results = Vector{Any}(undef, size(chosen_start, 2))

    count = 0
    for i in 1:size(chosen_start, 2)
        count += 1
        td = TwiceDifferentiable(obj, chosen_start[:, i], autodiff=:forward)
        res = optimize(
            td,
            chosen_start[:, i],
            #=Newton(),=#
            NewtonTrustRegion(;
                initial_delta=1e2,   # Δ₀
                delta_hat=1e4,  # max Δₖ
                eta=0.01,    # accept step when ρₖ > η
                rho_lower=0.1,   # shrink when ρₖ < ρ_lower
                rho_upper=0.9,   # grow   when ρₖ > ρ_upper
            ),
            Optim.Options(iterations=iters, f_abstol=tol, f_reltol=tol, g_abstol=NaN),
        )
        potential_results[i] = res
        if res.g_residual < 1.0
            break
        end
    end
    cut_results = potential_results[1:count]
    min_idx = argmin(res.minimum for res in cut_results)
    res = potential_results[min_idx]

    return (; res, td, num_iters, problem_starts)
end

function comovement_reg(data, dimvals, ranks; iters=300, tol=1e-05, num_starts=20, num_selected=5, spread=1)

    perm_mat = both_perm_mat(dimvals, ranks)
    perm_resp = (perm_mat*data)[:, 2:end]
    pred = data[:, 1:(end-1)]
    perm_cen = perm_resp .- mean(perm_resp, dims=2)

    res, td, num_iters, problem_starts = main_algorithm(perm_cen, pred, dimvals, ranks; iters, tol, num_starts, num_selected, spread)

    hess_non = hessian!(td, res.minimizer)
    hess_est = 0.5 .* (hess_non + hess_non')
    hess_eigs = real.(eigvals(hess_est))
    neg_eigs = hess_eigs[hess_eigs.<0.0]
    if !isempty(neg_eigs)
        @warn "Hessian has negative eigenvalues: $(neg_eigs). Using nearest positive semidefinite matrix."
        hess_pd = nearest_posdef(hess_est)
    else
        hess_pd = copy(hess_est)
    end
    theta_est = Optim.minimizer(res)

    delta_est, gamma_est, u3_est, u4_est, sigma_est =
        b_unpack_params(theta_est, dimvals, ranks)

    num_delta = ranks[1] * (dimvals[1] - ranks[1])
    stderrs = sqrt.(abs.(diag(inv(hess_pd))))
    delta_stderr = stderrs[1:num_delta]
    num_gamma = ranks[2] * (dimvals[2] - ranks[2])
    gamma_stderr = stderrs[(num_delta+1):(num_delta+num_gamma)]

    delta_star = delta_est[(dimvals[1]-ranks[1]+1):end, :]
    gamma_star = gamma_est[(dimvals[2]-ranks[2]+1):end, :]

    omega = omega_from_both(delta_star, gamma_star, dimvals, ranks)
    pi_mat = pi_from_both(u3_est, u4_est, dimvals, ranks)

    return (;
        res,
        delta_est,
        delta_stderr,
        gamma_est,
        gamma_stderr,
        u3_est,
        u4_est,
        sigma_est,
        hess_est,
        omega,
        pi_mat,
        num_iters,
        problem_starts,
    )
end

