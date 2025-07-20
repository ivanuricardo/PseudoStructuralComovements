
function unpack_params(theta, dimvals, ranks; p=1)
    num_delta = ranks[1] * (dimvals[1] - ranks[1])
    num_gamma = ranks[2] * (dimvals[2] - ranks[2])
    num_u3 = p * (ranks[1] * dimvals[1] - 1)
    num_u4 = p * (ranks[2] * dimvals[2])
    num_ll1 = Int(dimvals[1] * (dimvals[1] + 1) / 2) - 1
    num_ll2 = Int(dimvals[2] * (dimvals[2] + 1) / 2)

    @assert length(theta) == (num_delta + num_gamma + num_u3 + num_u4 + num_ll1 + num_ll2) "Parameter vector has wrong length!"

    idx_gamma = num_delta + num_gamma
    idx_u3 = num_delta + num_gamma + num_u3
    idx_u4 = num_delta + num_gamma + num_u3 + num_u4
    idx_ll1 = num_delta + num_gamma + num_u3 + num_u4 + num_ll1

    delta_vec = theta[1:num_delta]
    gamma_vec = theta[(num_delta+1):idx_gamma]
    u3_vec = theta[(idx_gamma+1):idx_u3]
    u4_vec = theta[(idx_u3+1):idx_u4]
    ll1_vec = theta[(idx_u4+1):idx_ll1]
    ll2_vec = theta[(idx_ll1+1):end]

    delta_star = reshape(delta_vec, ranks[1], dimvals[1] - ranks[1])
    delta = vcat(I, delta_star)

    gamma_star = reshape(gamma_vec, ranks[2], dimvals[2] - ranks[2])
    gamma = vcat(I, gamma_star)

    insert!(ll1_vec, 1, 1)
    ll1 = vec_to_ll(ll1_vec, dimvals[1])
    ll2 = vec_to_ll(ll2_vec, dimvals[2])

    if p == 1
        insert!(u3_vec, 1, 1)
        u3 = reshape(u3_vec, dimvals[1], ranks[1])
        u4 = reshape(u4_vec, dimvals[2], ranks[2])

        return (; delta, gamma, u3, u4, ll1, ll2)
    end

    t_u3 = eltype(u3_vec)
    t_u4 = eltype(u4_vec)
    u3 = zeros(t_u3, p * dimvals[1], ranks[1])
    u4 = zeros(t_u4, p * dimvals[2], ranks[2])
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
    return (; delta, gamma, u3, u4, ll1, ll2)

end

function pack_params(delta_star, gamma_star, u3, u4, ll1, ll2; p=1)
    n1, n2 = size(u3)
    true_n1 = n1 ÷ p
    vec_u3 = vecb(u3, true_n1)
    removek!(vec_u3, true_n1 * n2 - 1)
    vec_ll1 = vech(ll1)
    removek!(vec_ll1, 1)
    vec_ll2 = vech(ll2)

    return vcat(vec(delta_star), vec(gamma_star), vec_u3, vec(u4), vec_ll1, vec_ll2)
end

function rand_init(dimvals, ranks; p=1)
    n1_r1 = dimvals[1] - ranks[1]
    n2_r2 = dimvals[2] - ranks[2]

    coef = generate_rrmar_coef(dimvals, ranks; p, maxeigen=0.9)

    delta = coef.delta
    gamma = coef.gamma
    delta_init = delta[(n1_r1+1):end, :]
    gamma_init = gamma[(n2_r2+1):end, :]
    return pack_params(delta_init, gamma_init, coef.u3, coef.u4, I(dimvals[1]), I(dimvals[2]); p)

end

function init_alg(resp, pred, dimvals, ranks; p=1)

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

    perm_mat = perm_matrix(dimvals, ranks)
    omega = create_omega(delta_star, gamma_star, dimvals, ranks)
    pi_mat = omega * perm_mat * coef
    pi_star = pi_mat[(pdims-r+1):end, :]
    u4_est, u3_est = nearest_kron(pi_star', size(u2), size(u1))

    s = u3_est[1, 1]

    u3 = u3_est / s
    u4 = u4_est * s

    if p != 1
        u3_bot = randn(dimvals[1], ranks[1])
        s_alt = u3_bot[1, 1]
        u3_bot = copy(u3_bot) / s_alt
        u4_bot = randn(dimvals[2], ranks[2]) * s_alt
        u3 = vcat(u3, u3_bot)
        u4 = vcat(u4, u4_bot)
    end

    return pack_params(delta_star, gamma_star, u3, u4, I(dimvals[1]), I(dimvals[2]); p)

end

function loglike(theta, resp, pred, dimvals, ranks; p=1)
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    delta_rot, gamma_rot, u3, u4, ll1, ll2 = unpack_params(theta, dimvals, ranks; p)
    sigma1 = ll1 * ll1'
    sigma2 = ll2 * ll2'
    ll = kron(ll2, ll1)

    delta_star = delta_rot[(N1_r1+1):end, :]
    gamma_star = gamma_rot[(N2_r2+1):end, :]

    pi_mat = create_pi(u3, u4, dimvals, ranks; p)

    obs = size(resp, 2)
    omega = create_omega(delta_star, gamma_star, dimvals, ranks)
    if p > 1
        omega_tilde, pi_tilde, ll = make_companion(omega, pi_mat; ll)
        sparse_omega = sparse(omega_tilde)
        sparse_pi = sparse(pi_tilde)
    else
        sparse_omega = sparse(omega)
        sparse_pi = sparse(pi_mat)
    end

    X = sparse_omega * ll
    precision_matrix = inv(X') * inv(X)
    sparse_precision = sparse(precision_matrix)

    sse = 0.0

    for i = 2:obs
        yt = @view resp[:, i]
        yt_m1 = @view pred[:, i]
        resid = sparse_omega * yt - sparse_pi * yt_m1
        sse += dot(resid, sparse_precision * resid)
    end

    det_term1 = det(sigma1)
    det_term2 = det(sigma2)
    if det_term1 <= 0.0 || det_term2 <= 0.0
        # quadratic penalty on negative/zero dets
        pen = 1e6 * ((min(det_term1, 0.0))^2 + (min(det_term2, 0.0))^2)
        return 0.5 * sse + pen
    end

    # no need for omegas because det = 1
    logdet_term1 = dimvals[2] * log(det_term1)
    logdet_term2 = dimvals[1] * log(det_term2)

    return 0.5 * ((obs - 1) * (logdet_term1 + logdet_term2) + sse)
end

function comovement_init(resp, pred, dimvals, ranks; iters=5, tol=1e-08, num_starts=50, num_selected=10, p=1)
    some_init = init_alg(resp, pred, dimvals, ranks; p)
    init_length = length(some_init)
    potential_starts = fill(NaN, init_length + 1, num_starts)
    obj = tet -> loglike(tet, resp, pred, dimvals, ranks; p)

    for i in 1:num_starts
        if i == 1
            both_init = copy(some_init)
        else
            both_init = rand_init(dimvals, ranks; p)
        end
        potential_starts[1:(end-1), i] = both_init
        td = TwiceDifferentiable(obj, both_init, autodiff=:forward)

        res = optimize(
            td,
            both_init,
            BFGS(),
            Optim.Options(iterations=iters, f_abstol=tol, f_reltol=tol),
        )
        potential_starts[end, i] = res.minimum
        potential_starts[1:(end-1), i] = res.minimizer
    end
    chosen_idx = partialsortperm(potential_starts[end, :], 1:num_selected)
    chosen_start = potential_starts[1:(end-1), chosen_idx]

    return chosen_start

end

function check_neg_eigs(td, res)
    hess_non = hessian!(td, res.minimizer)
    hess_est = 0.5 .* (hess_non + hess_non')
    hess_eigs = real.(eigvals(hess_est))
    neg_eigs = hess_eigs[hess_eigs.<0.0]
    if any(neg_eigs .< 1e-01)
        return true
    end
    return false
end

function main_algorithm(resp, pred, dimvals, ranks; iters=1000, tol=1e-08, num_starts=50, num_selected=10, p=1, grad_tol=1e-01)
    obj = tet -> loglike(tet, resp, pred, dimvals, ranks; p)
    chosen_start = comovement_init(resp, pred, dimvals, ranks; iters=5, tol=grad_tol, num_starts, num_selected, p)
    potential_results = []
    td = nothing
    res = nothing

    # Run all optimizations and collect results
    for i in 1:size(chosen_start, 2)
        td = TwiceDifferentiable(obj, chosen_start[:, i], autodiff=:forward)
        res = optimize(
            td,
            chosen_start[:, i],
            BFGS(),
            Optim.Options(iterations=iters, f_abstol=tol, f_reltol=tol, g_abstol=grad_tol),
        )
        neg_eig_check = check_neg_eigs(td, res)
        push!(potential_results, (res, td, neg_eig_check))
    end

    # Filter results with low gradient norm
    valid_results = [
        r for (r, td, neg_check) in potential_results
        if !neg_check && r.g_residual < 1.0
    ]

    if !isempty(valid_results)
        # Choose valid result with smallest objective value
        min_obj_idx = argmin([r.minimum for r in valid_results])
        res = valid_results[min_obj_idx]
    else
        # Fallback: choose result with smallest gradient norm value across all runs
        min_obj_idx = argmin([r[1].g_residual for r in potential_results])
        res, td = potential_results[min_obj_idx][1:2]  # Extract res and td
    end

    return (; res, td)
end

function comovement_reg(data, dimvals, ranks; iters=1000, tol=1e-08, num_starts=50, num_selected=10, p=1)

    if p != 1
        if prod(dimvals) * p != size(data, 1)
            data = companion_data(data, p)
        end
        perm_mat = kron(I(p), perm_matrix(dimvals, ranks))
    else
        perm_mat = perm_matrix(dimvals, ranks)
    end
    perm_resp = (perm_mat*data)[:, 2:end]
    pred = data[:, 1:(end-1)]
    resp = perm_resp .- mean(perm_resp, dims=2)

    all_results = []
    count = 0

    # First optimization attempt
    res, td = main_algorithm(resp, pred, dimvals, ranks; iters, tol, num_starts, num_selected, p)
    first_neg_eig_check = check_neg_eigs(td, res)
    push!(all_results, (res, td, first_neg_eig_check))

    # Subsequent attempts if needed
    while res.g_residual > 1e-01 || first_neg_eig_check
        res, td = main_algorithm(resp, pred, dimvals, ranks; iters, tol, num_starts, num_selected, p)
        next_neg_eig_check = check_neg_eigs(td, res)
        push!(all_results, (res, td, next_neg_eig_check))
        count += 1
        (count >= 5) && break  # Max 5 additional attempts
    end

    # Select best result: prioritize valid (no neg eigs + low grad) then fallback to min objective
    valid_results = [
        r for (r, td, neg_check) in all_results
        if !neg_check && r.g_residual < 1.0
    ]

    if !isempty(valid_results)
        # Choose valid result with smallest objective value
        min_obj_idx = argmin([r.minimum for r in valid_results])
        res = valid_results[min_obj_idx]
        println("Selected valid result (g_res=$(res.g_residual))")
    else
        # Fallback: choose result with smallest objective value across all runs
        min_obj_idx = argmin([r[1].minimum for r in all_results])
        res, td = all_results[min_obj_idx][1:2]  # Extract res and td
        @warn "No valid results! Using fallback (obj=$(res.minimum), g_res=$(res.g_residual))"
    end

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

    delta_est, gamma_est, u3_est, u4_est, ll1_est, ll2_est =
        unpack_params(theta_est, dimvals, ranks; p)

    num_delta = ranks[1] * (dimvals[1] - ranks[1])
    stderrs = sqrt.(abs.(diag(inv(hess_pd))))
    delta_stderr = stderrs[1:num_delta]
    num_gamma = ranks[2] * (dimvals[2] - ranks[2])
    gamma_stderr = stderrs[(num_delta+1):(num_delta+num_gamma)]

    delta_star = delta_est[(dimvals[1]-ranks[1]+1):end, :]
    gamma_star = gamma_est[(dimvals[2]-ranks[2]+1):end, :]

    omega = create_omega(delta_star, gamma_star, dimvals, ranks)
    pi_mat = create_pi(u3_est, u4_est, dimvals, ranks; p)
    sigma1_est = ll1_est * ll1_est'
    sigma2_est = ll2_est * ll2_est'

    return (;
        res,
        delta_est,
        delta_stderr,
        gamma_est,
        gamma_stderr,
        u3_est,
        u4_est,
        sigma1_est,
        sigma2_est,
        hess_est,
        omega,
        pi_mat,
    )
end

