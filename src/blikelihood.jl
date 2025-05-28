
function b_unpack_params(theta, dimvals, ranks; p=1)
    pdims = prod(dimvals)
    num_delta = ranks[1] * (dimvals[1] - ranks[1])
    num_gamma = ranks[2] * (dimvals[2] - ranks[2])
    num_u3 = p * (ranks[1] * dimvals[1] - 1)
    num_u4 = p * (ranks[2] * dimvals[2])
    num_sigma = Int(pdims * (pdims + 1) / 2)

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
    delta = vcat(I, delta_star)

    gamma_star = reshape(gamma_vec, ranks[2], dimvals[2] - ranks[2])
    gamma = vcat(I, gamma_star)
    ll = vec_to_ll(ll_vec, pdims)

    if p == 1
        insert!(u3_vec, 1, 1)
        u3 = reshape(u3_vec, dimvals[1], ranks[1])
        u4 = reshape(u4_vec, dimvals[2], ranks[2])

        return (; delta, gamma, u3, u4, ll)
    end

    u3 = zeros(p * dimvals[1], ranks[1])
    u4 = zeros(p * dimvals[2], ranks[2])
    insertk!(u3_vec, num_u3 ÷ p)
    for k in 0:(p-1)

        i1, i2 = k * dimvals[1] + 1, (k + 1) * dimvals[1]
        j1, j2 = k * dimvals[2] + 1, (k + 1) * dimvals[2]

        u3_prod = dimvals[1] * ranks[1]
        u4_prod = dimvals[2] * ranks[2]
        k1, k2 = k * u3_prod + 1, (k + 1) * u3_prod
        l1, l2 = k * u4_prod + 1, (k + 1) * u4_prod

        u3[i1:i2, :] .= reshape(u3_vec[k1:k2], dimvals[1], ranks[1])
        u4[j1:j2, :] .= reshape(u4_vec[l1:l2], dimvals[2], ranks[2])
    end
    return (; delta, gamma, u3, u4, ll)

end

function b_pack_params(delta, gamma, u3, u4, ll; p=1)
    n1, n2 = size(u3)
    true_n1 = n1 ÷ p
    vec_u3 = vecb(u3, true_n1)
    removek!(vec_u3, true_n1 * n2 - 1)
    return vcat(vec(delta), vec(gamma), vec_u3, vec(u4), vech(ll))
end

function rand_init(dimvals, ranks; p=1)
    n1_r1 = dimvals[1] - ranks[1]
    n2_r2 = dimvals[2] - ranks[2]

    coef = generate_rrmar_coef(dimvals, ranks; maxeigen=0.9)

    delta = coef.delta
    gamma = coef.gamma
    delta_init = delta[(n1_r1+1):end, :]
    gamma_init = gamma[(n2_r2+1):end, :]
    return b_pack_params(delta_init, gamma_init, coef.u3, coef.u4, I(prod(dimvals)); p)

end

function init_both(resp, pred, dimvals, ranks; pack_params=true, p=1)

    pdims = prod(dimvals)
    r = prod(ranks)
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    coef = ols_coef(resp, pred)
    if p != 1
        coef = coef[1:pdims, 1:pdims]
    end

    tensor_coef =
        matten(coef, [1, 2], [3, 4], [dimvals[1], dimvals[2], dimvals[1], dimvals[2]])
    flat1 = tenmat(tensor_coef, row=1)
    flat2 = tenmat(tensor_coef, row=2)
    u1 = svd(flat1).U[:, 1:ranks[1]]
    u2 = svd(flat2).U[:, 1:ranks[2]]

    delta = nullspace(u1')
    rotate_u!(delta)
    delta_star = delta[(N1_r1+1):end, :]
    gamma = nullspace(u2')
    rotate_u!(gamma)
    gamma_star = gamma[(N2_r2+1):end, :]

    perm_mat = both_perm_mat(dimvals, ranks)
    omega = omega_from_both(delta_star, gamma_star, dimvals, ranks)
    pi_mat = omega * perm_mat * coef
    pi_star = pi_mat[(pdims-r+1):end, :]
    u4_est, u3_est = nearest_kron(pi_star', size(u2), size(u1))

    s = u3_est[1, 1]

    u3 = u3_est / s
    u4 = u4_est * s

    if p != 1
        u3_bot = randn(dimvals[1], ranks[1])
        u4_bot = randn(dimvals[2], ranks[2])
        u3 = vcat(u3, u3_bot)
        u4 = vcat(u4, u4_bot)
    end

    if pack_params
        return b_pack_params(delta_star, gamma_star, u3, u4, I(prod(dimvals)); p)
    else
        return (; u1, u2, u3, u4)
    end

end

function both_loglike(theta, resp, pred, dimvals, ranks; p=1)
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    delta_rot, gamma_rot, u3, u4, ll = b_unpack_params(theta, dimvals, ranks; p)
    delta_star = delta_rot[(N1_r1+1):end, :]
    gamma_star = gamma_rot[(N2_r2+1):end, :]

    pi_mat = pi_from_both(u3, u4, dimvals, ranks; p)

    obs = size(resp, 2)
    omega = omega_from_both(delta_star, gamma_star, dimvals, ranks)
    if p > 1
        omega_tilde, pi_tilde, ll = make_companion(omega, pi_mat, ll)
        sparse_omega = sparse(omega_tilde)
        sparse_pi = sparse(pi_tilde)
    else
        sparse_omega = sparse(omega)
        sparse_pi = sparse(pi_mat)
    end
    det_term = det(0.5 .* (ll * ll' + ll * ll'))
    if det_term <= 0.0
        return 1e10
    end

    logdet_term = log(det_term)
    X = sparse_omega * ll
    precision_matrix = inv(X') * inv(X)

    sse = 0.0

    for i = 2:obs
        yt = @view resp[:, i]
        yt_m1 = @view pred[:, i]
        resid = sparse_omega * yt - sparse_pi * yt_m1
        sse += dot(resid, precision_matrix * resid)
    end

    return 0.5 * ((obs - 1) * logdet_term + sse)
end

function comovement_init(resp, pred, dimvals, ranks; iters=5, tol=1e-05, num_starts=20, num_selected=10)
    some_init = init_both(resp, pred, dimvals, ranks)
    init_length = length(some_init)
    potential_starts = fill(NaN, init_length + 1, num_starts)
    obj = tet -> both_loglike(tet, resp, pred, dimvals, ranks)

    for i in 1:num_starts
        if i == 1
            both_init = copy(some_init)
        else
            #=both_init = copy(some_init) .+ rand_init(dimvals, ranks)=#
            both_init = rand_init(dimvals, ranks)
        end
        potential_starts[1:(end-1), i] = both_init
        td = TwiceDifferentiable(obj, both_init, autodiff=:forward)

        res = optimize(
            td,
            both_init,
            #=LBFGS(),=#
            NewtonTrustRegion(;
                initial_delta=1e3,      # start near the 2.2e3 step that first reduced f
                delta_hat=5e4,      # cap region to ~5×10^3 so you don’t waste work on huge proposals
                eta=0.02,     # accept any step with ρₖ > 0.05 (instead of 0.1)
                rho_lower=0.05,     # only shrink if ρₖ < 0.05 (more forgiving)
                rho_upper=0.5,      # grow if ρₖ > 0.7 (easier to expand when things look good)
            ),
            Optim.Options(iterations=iters, f_abstol=tol, f_reltol=tol),
        )
        potential_starts[end, i] = res.minimum
        potential_starts[1:(end-1), i] = res.minimizer
    end
    chosen_idx = partialsortperm(potential_starts[end, :], 1:num_selected)
    chosen_start = potential_starts[1:(end-1), chosen_idx]

    return chosen_start

end

function main_algorithm(resp, pred, dimvals, ranks; iters=100, tol=1e-05, num_starts=20, num_selected=10)

    obj = tet -> both_loglike(tet, resp, pred, dimvals, ranks)
    td = nothing
    res = nothing

    chosen_start = comovement_init(resp, pred, dimvals, ranks; iters=5, tol=1e-01, num_starts, num_selected)
    potential_results = Vector{Any}(undef, size(chosen_start, 2))

    count = 0
    for i in 1:size(chosen_start, 2)
        count += 1
        td = TwiceDifferentiable(obj, chosen_start[:, i], autodiff=:forward)
        res = optimize(
            td,
            chosen_start[:, i],
            #=LBFGS(),=#
            NewtonTrustRegion(
                initial_delta=1e3,      # start near the 2.2e3 step that first reduced f
                delta_hat=5e4,      # cap region to ~5×10^3 so you don’t waste work on huge proposals
                eta=0.02,     # accept any step with ρₖ > 0.05 (instead of 0.1)
                rho_lower=0.05,     # only shrink if ρₖ < 0.05 (more forgiving)
                rho_upper=0.5,      # grow if ρₖ > 0.7 (easier to expand when things look good)
            ),
            Optim.Options(iterations=iters, f_abstol=tol, f_reltol=tol, g_abstol=1e-01, store_trace=true),
        )
        potential_results[i] = res
        if res.g_residual < 1e-01
            break
        end
    end
    cut_results = potential_results[1:count]
    min_idx = argmin(res.minimum for res in cut_results)
    res = potential_results[min_idx]
    println(res.iterations)

    return (; res, td)
end

function comovement_reg(data, dimvals, ranks; iters=100, tol=1e-05, num_starts=20, num_selected=10, p=1)

    if p != 1
        data = companion_data(data, p)
        perm_mat = kron(I(p), both_perm_mat(dimvals, ranks))
    else
        perm_mat = both_perm_mat(dimvals, ranks)
    end
    perm_resp = (perm_mat*data)[:, 2:end]
    pred = data[:, 1:(end-1)]
    resp = perm_resp .- mean(perm_resp, dims=2)

    res, td = main_algorithm(resp, pred, dimvals, ranks; iters, tol, num_starts, num_selected)

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
    )
end

