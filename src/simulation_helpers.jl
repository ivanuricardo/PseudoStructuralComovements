
function round!(A; digits=0)
    A .= round.(A; digits=digits)
end

function isstable(A, maxeigen)
    stab_cond = maximum(abs.(eigvals(A)))
    return stab_cond < maxeigen
end

function generate_rrmar_coef(dimvals, ranks; p=1, maxeigen=0.9, coef_scale = 1)

    A = fill(NaN, p * prod(dimvals), p * prod(dimvals))
    delta = fill(NaN, dimvals[1], dimvals[1] - ranks[1])
    gamma = fill(NaN, dimvals[2], dimvals[2] - ranks[2])
    u3 = fill(NaN, p * dimvals[1], ranks[1])
    u4 = fill(NaN, p * dimvals[2], ranks[2])
    stabit = 0

    while true
        stabit += 1
        delta_star = randn(ranks[1], dimvals[1] - ranks[1])
        delta .= vcat(I, delta_star)
        gamma_star = randn(ranks[2], dimvals[2] - ranks[2])
        gamma .= vcat(I, gamma_star)
        omega_star = create_omega(delta_star, gamma_star, dimvals, ranks)
        u3_scale = fill(NaN, p)
        count = 0

        for i in 1:dimvals[1]:(p*dimvals[1])
            count += 1
            u3_range = i:i+dimvals[1]-1
            u3_partial = randn(dimvals[1], ranks[1])
            u3_scale[count] = u3_partial[1, 1]
            u3[u3_range, :] .= u3_partial / u3_scale[count]
        end
        count = 0
        for i in 1:dimvals[2]:(p*dimvals[2])
            count += 1
            u4_range = i:i+dimvals[2]-1
            u4_partial = coef_scale .* randn(dimvals[2], ranks[2])
            u4[u4_range, :] = u4_partial * u3_scale[count]
        end
        pi_mat = create_pi(u3, u4, dimvals, ranks; p)
        perm_mat = perm_matrix(dimvals, ranks)
        if p == 1
            # given Ω^* P y_t = Π y_{t-1} + Ω^* P e_t
            # A^* = inv(Ω^*) Π = inv(Ω P') Π 
            # A = P' Ω^* Π to obtain y_t = A y_{t-1} + e_t
            A .= perm_mat' * inv(omega_star) * pi_mat
        else
            omega_tilde, pi_tilde = make_companion(omega_star, pi_mat)
            large_perm = kron(I(p), perm_mat)
            A .= large_perm' * inv(omega_tilde) * pi_tilde
        end
        coef_eigs = round.(sort(abs.(eigvals(A)), rev=true), digits=6)
        last_idx = findlast(!iszero, coef_eigs)
        if isnothing(last_idx)
            @error "No stable coefficients found!"
        end

        if isstable(A, maxeigen) && (coef_eigs[last_idx] > 0.4)
            break
        elseif stabit > 1e4
            coef_scale *= 0.1
            stabit = 0
        end
    end

    sorted_eigs = sort(abs.(eigvals(A)), rev=true)
    return (; A, delta, gamma, u3, u4, stabit, sorted_eigs, coef_scale)
end

function simulate_rrmar_data(
    dimvals::AbstractVector,
    ranks::AbstractVector,
    obs::Int;
    A=nothing,
    snr=0.7,
    burnin=100,
    matrix_err=false,
    p=1
)
    if isnothing(A)
        A = generate_rrmar_coef(dimvals, ranks; p)
    end

    dim1, dim2 = dimvals[1], dimvals[2]
    coef = A.A
    rho = maximum(abs.(eigvals(coef)))
    sigma_val = rho / snr
    diagerr = repeat([sigma_val], dim1 * dim2)
    sigma = diagm(diagerr)
    pre_sigma2, pre_sigma1 = nearest_kron(sigma, (dim2, dim2), (dim1, dim1))
    sigma1 = abs.(pre_sigma1) / abs(pre_sigma1[1])
    sigma2 = abs.(pre_sigma2) * abs(pre_sigma1[1])
    pre_data = zeros(p * dimvals[1] * dimvals[2], obs)
    round!(sigma1; digits=15)
    round!(sigma2; digits=15)

    if matrix_err
        d = MatrixNormal(zeros(dim1, dim2), sigma1, sigma2)
        mat_errs = rand(d, obs)
        vec_errs = hcat([vec(mat_errs[i]) for i in eachindex(mat_errs)]...)

        for i = 2:obs
            pre_data[:, i] .= coef * pre_data[:, i-1] + vec_errs[:, i]
        end
        data = pre_data[:, (burnin+1):end]
        return (; data, coef, sigma1, sigma2)
    else
        d = MultivariateNormal(zeros(dimvals[1] * dimvals[2]), diagerr)
        if p == 1
            for i = 2:obs
                vec_epsilon = rand(d)
                pre_data[:, i] .= coef * pre_data[:, i-1] + vec_epsilon
            end
            data = pre_data[:, (burnin+1):end]
            sigma = copy(diagerr)
            return (; data, coef, sigma)
        end

        for i = 2:obs
            vec_epsilon = vcat(rand(d), zeros(prod(dimvals) * (p - 1)))
            pre_data[:, i] .= coef * pre_data[:, i-1] + vec_epsilon
        end
        data = pre_data[:, (burnin+1):end]
        sigma = copy(diagerr)
        return (; data, coef, sigma)
    end
end
