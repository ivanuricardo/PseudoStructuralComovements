function perturb!(A; digs=10)
    @inbounds for i in eachindex(A)
        A[i] += 10^(-float(digs)) * randn()
    end
    return A
end

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
    vec_u3 = filter(x -> !isapprox(x, 1.0; atol=1e-12), vec(u3))
    vec_ll1 = vech(ll1)
    deleteat!(vec_ll1, 1)
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

function stack_with_random(A::AbstractMatrix{T}, p::Integer) where T<:Real
    N, r = size(A)
    mats = Vector{Matrix{T}}(undef, p+1)
    mats[1] = A
    for i in 1:p
        B = randn(T, N, r)
        B[1,1] = A[1,1]       # enforce the (1,1)–entry
        mats[i+1] = B
    end
    return vcat(mats...)      # stack them all vertically
end

function init_alg(data, dimvals, ranks; p=1)

    data = data[1:prod(dimvals), :]
    obs = size(data, 2)
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]
    ten_data = reshape(data, (dimvals[1], dimvals[2], obs))
    perm_data = permutedims(ten_data, (3, 1, 2))
    (u1_est, u2_est, u3_est, u4_est, sigma1, sigma2, s1, s2) = lock(R_LOCK) do
        rrmar_init = R"""
        d1 = $(dimvals[1])
        d2 = $(dimvals[2])
        r1 = $(ranks[1])
        r2 = $(ranks[2])
        data <- $(perm_data)
        est <- matAR.RR.est(data, method = "RRMLE", k1 = r1, k2 = r2, tol = 1e-08)
        s1 <- t(est$loading$U1) %*% est$A1 %*% est$loading$V1
        s2 <- t(est$loading$U2) %*% est$A2 %*% est$loading$V2
        """
        @rget est s1 s2
        (
            est[:loading][:U1][:, 1:ranks[1]],   # u1_est
            est[:loading][:U2][:, 1:ranks[2]],   # u2_est
            est[:loading][:V1][:, 1:ranks[1]],   # u3_est
            est[:loading][:V2][:, 1:ranks[2]],   # u4_est
            est[:Sig1],                           # sigma1
            est[:Sig2],                            # sigma2
            s1,
            s2
        )
    end

    q1 = inv(u1_est[(N1_r1+1):end, :])
    q2 = inv(u2_est[(N2_r2+1):end, :])
    u1_star = u1_est * q1
    u2_star = u2_est * q2
    u3_star = u3_est * s1[1:ranks[1], 1:ranks[1]] * inv(q1)'
    u4_star = u4_est * s2[1:ranks[2], 1:ranks[2]] * inv(q2)'

    delta_star = -u1_star[1:N1_r1, :]'
    gamma_star = -u2_star[1:N2_r2, :]'

    u3 = u3_star / u3_star[1]
    u4 = u4_star * u3_star[1]

    chol1 = cholesky(Symmetric(sigma1))
    ll1 = chol1.L
    chol2 = cholesky(Symmetric(sigma2))
    ll2 = chol2.L

    ss = ll1[1, 1]
    ll1_rot = ll1 / ss
    ll2_rot = ll2 * ss

    if p != 1
        u3 = stack_with_random(u3, p-1)
        u4 = stack_with_random(u4, p-1)
    end

    return pack_params(delta_star, gamma_star, u3, u4, ll1_rot, ll2_rot; p)

end

function loglike_calc(theta, resp, pred, dimvals, ranks; p=1)
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    delta_rot, gamma_rot, u3, u4, ll1, ll2 = unpack_params(theta, dimvals, ranks; p)
    sigma1 = ll1 * ll1'
    sigma2 = ll2 * ll2'
    ll = kron(ll2, ll1)
    if p == 1
        perm_mat = perm_matrix(dimvals, ranks)
    else
        z = zeros(prod(dimvals), (p-1) * prod(dimvals))
        perm_mat = [ perm_matrix(dimvals, ranks) z; z' I((p-1) * prod(dimvals)) ]
    end

    delta_star = delta_rot[(N1_r1+1):end, :]
    gamma_star = gamma_rot[(N2_r2+1):end, :]

    pi_mat = create_pi(u3, u4, dimvals, ranks; p)

    obs = size(resp, 2)
    omega = create_omega(delta_star, gamma_star, dimvals, ranks)
    if p > 1
        omega_tilde, pi_tilde, ll = make_companion(omega, pi_mat; ll)
        sparse_omega = sparse(omega_tilde * perm_mat)
        sparse_pi = sparse(pi_tilde)
    else
        sparse_omega = sparse(omega * perm_mat)
        sparse_pi = sparse(pi_mat)
    end

    X = sparse_omega * ll
    Xfact = try
        lu(X)
    catch e
        isa(e, SingularException) || rethrow(e)
        return 1e9
    end

    resid_buf = zeros(eltype(theta), size(X, 1))
    z_buf     = similar(resid_buf)
    sse = 0
    for i in 1:obs
        yt   = @view resp[:, i]
        yt_m1 = @view pred[:, i]
        mul!(resid_buf, sparse_omega, yt)
        mul!(resid_buf, sparse_pi, yt_m1, -1.0, 1.0)
        ldiv!(z_buf, Xfact, resid_buf)
        sse += dot(z_buf, z_buf)
    end

    d1 = diag(ll1)
    d2 = diag(ll2)
    if any(x -> x <= 0.0, d1) || any(x -> x <= 0.0, d2)
        pen = 1e6 * (sum(x -> min(x, 0.0)^2, d1) + sum(x -> min(x, 0.0)^2, d2))
        return 0.5 * sse + pen
    end
    logdet_term1 = dimvals[2] * 2 * sum(log, d1)   # log(det(sigma)) = 2 * sum(log(diag(L)))
    logdet_term2 = dimvals[1] * 2 * sum(log, d2)

    return 0.5 * ((obs) * (logdet_term1 + logdet_term2) + sse)
end

function loglike(theta, resp, pred, dimvals, ranks; p=1)
    ell = try
        loglike_calc(theta, resp, pred, dimvals, ranks; p)
    catch e
        if isa(e, InterruptException) ||
           isa(e, OutOfMemoryError) ||
           isa(e, StackOverflowError) ||
           isa(e, MethodError) ||
           isa(e, UndefVarError) ||
           isa(e, LoadError) ||
           isa(e, InitError)
            rethrow(e)
        end
        @warn "Log Likelihood returns an error!" exception=(e, catch_backtrace())
        return 1e12
    end

    if !isfinite(ell) || abs(ell) > 1e12
        return 1e12
    end

    return ell
end

function comovement_init(data, resp, pred, dimvals, ranks; iters=5, tol=1e-10, num_starts=100, num_selected=10, p=1)
    some_init = init_alg(data, dimvals, ranks; p)
    init_length = length(some_init)
    potential_starts = fill(NaN, init_length + 1, num_starts)
    obj = tet -> loglike(tet, resp, pred, dimvals, ranks; p)

    for i in 1:num_starts
        both_init = copy(some_init) .+ 0.01 .* randn(length(some_init))
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
    if any(neg_eigs .< -1e-02)
        return true
    end
    return false
end

function main_algorithm(data, resp, pred, dimvals, ranks; iters=1000, tol=1e-10, num_starts=100, num_selected=10, p=1, grad_tol=1e-02)
    obj = tet -> loglike(tet, resp, pred, dimvals, ranks; p)
    chosen_start = comovement_init(data, resp, pred, dimvals, ranks; iters=5, tol=grad_tol, num_starts, num_selected, p)
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

function comovement_reg(data, dimvals, ranks; iters=1000, tol=1e-10, num_starts=100, num_selected=10, p=1)

    nrows = size(data, 1)
    expected = prod(dimvals) * p
    if p != 1 && expected != nrows
        data = companion_data(data, p)
    end

    perm_resp = data[:, 2:end]
    pred = data[:, 1:(end-1)] .- mean(data[:, 1:(end-1)])
    resp = perm_resp .- mean(perm_resp, dims=2)

    all_results = []
    count = 0

    # First optimization attempt
    res, td = main_algorithm(data, resp, pred, dimvals, ranks; iters, tol, num_starts, num_selected, p)
    first_neg_eig_check = check_neg_eigs(td, res)
    push!(all_results, (res, td, first_neg_eig_check))

    while res.g_residual > 1e-02 || first_neg_eig_check
        res, td = main_algorithm(data, resp, pred, dimvals, ranks; iters, tol, num_starts, num_selected, p)
        next_neg_eig_check = check_neg_eigs(td, res)
        push!(all_results, (res, td, next_neg_eig_check))
        count += 1
        (count >= 0) && break  # Max 3 additional attempts
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
    else
        # Fallback: choose result with smallest objective value across all runs
        min_obj_idx = argmin([r[1].minimum for r in all_results])
        res, td = all_results[min_obj_idx][1:2]  # Extract res and td
        theta_est = Optim.minimizer(res)

        delta_est, gamma_est, u3_est, u4_est, ll1_est, ll2_est = unpack_params(theta_est, dimvals, ranks; p)
        @warn "No valid results! Using fallback, g_res=$(res.g_residual), delta is $(delta_est), gamma is $(gamma_est)"
    end

    hess_non = hessian!(td, res.minimizer)
    hess_est = 0.5 .* (hess_non + hess_non')
    hess_eigs = real.(eigvals(hess_est))

    max_eig = maximum(abs.(hess_eigs))
    min_eig = minimum(abs.(hess_eigs))

    neg_eigs = hess_eigs[hess_eigs.<0.0]
    if !isempty(neg_eigs)
        @warn "Hessian has negative eigenvalues: $(neg_eigs). Using nearest positive semidefinite matrix."
        hess_pd = nearest_posdef(hess_est)
    elseif min_eig < 1e-09
        @warn "Hessian is near singular! min eig is $min_eig and condition number is $(cond(hess_est))"
        hess_pd = nearest_posdef(hess_est)
    else
        hess_pd = copy(hess_est)
    end

    num_delta = ranks[1] * (dimvals[1] - ranks[1])
    stderrs = sqrt.(abs.(diag(inv(hess_pd))))
    delta_stderr = stderrs[1:num_delta]
    num_gamma = ranks[2] * (dimvals[2] - ranks[2])
    gamma_stderr = stderrs[(num_delta+1):(num_delta+num_gamma)]

    theta_est = Optim.minimizer(res)

    delta_est, gamma_est, u3_est, u4_est, ll1_est, ll2_est = unpack_params(theta_est, dimvals, ranks; p)

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

function rrmar(data, dimvals, ranks)

    resp_before_demean = data[:, 2:end]
    pred_before_demean = data[:, 1:end-1]
    resp = resp_before_demean .- mean(resp_before_demean, dims = 2)
    pred = pred_before_demean .- mean(pred_before_demean, dims = 2)
    theta_est = init_alg(data, dimvals, ranks)
    delta_est, gamma_est, u3_est, u4_est, ll1_est, ll2_est = unpack_params(theta_est, dimvals, ranks; p=1)
    ll = loglike(theta_est, resp, pred, dimvals, ranks)

    delta_star = delta_est[(dimvals[1]-ranks[1]+1):end, :]
    gamma_star = gamma_est[(dimvals[2]-ranks[2]+1):end, :]

    omega = create_omega(delta_star, gamma_star, dimvals, ranks)
    pi_mat = create_pi(u3_est, u4_est, dimvals, ranks; p=1)
    sigma1_est = ll1_est * ll1_est'
    sigma2_est = ll2_est * ll2_est'

    return (;
        delta_est,
        gamma_est,
        u3_est,
        u4_est,
        sigma1_est,
        sigma2_est,
        omega,
        pi_mat,
        ll,
    )

end
