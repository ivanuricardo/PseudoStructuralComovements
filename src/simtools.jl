
function round!(A; digits=0)
    A .= round.(A; digits=digits)
end

function isstable(A, maxeigen)
    stab_cond = maximum(abs.(eigvals(A)))
    return stab_cond < maxeigen
end

function generate_rrmar_coef(dimvals, ranks; p=1, maxeigen=0.9)

    A = fill(NaN, prod(dimvals), prod(dimvals))
    delta = fill(NaN, dimvals[1], dimvals[1] - ranks[1])
    gamma = fill(NaN, dimvals[2], dimvals[2] - ranks[2])
    u3 = fill(NaN, p * dimvals[1], ranks[1])
    u4 = fill(NaN, p * dimvals[2], ranks[2])
    stabit = 0

    u3_cont = Vector{Any}(undef, p)
    u4_cont = Vector{Any}(undef, p)
    while true
        stabit += 1
        delta_star = randn(ranks[1], dimvals[1] - ranks[1])
        delta .= vcat(I, delta_star)
        gamma_star = randn(ranks[2], dimvals[2] - ranks[2])
        gamma .= vcat(I, gamma_star)
        omega = omega_from_both(delta_star, gamma_star, dimvals, ranks)

        u3 = randn(dimvals[1], ranks[1])
        u3_scale = u3[1, 1]
        u3 = copy(u3) / u3_scale
        u4 = randn(dimvals[2], ranks[2])
        u4 = copy(u4) * u3_scale
        pi_mat = pi_from_both(u3, u4, dimvals, ranks)
        perm_mat = both_perm_mat(dimvals, ranks)
        A .= inv(omega * perm_mat) * pi_mat

        if isstable(A, maxeigen)
            break
        end
    end

    sorted_eigs = sort(abs.(eigvals(A)), rev=true)
    return (; A, delta, gamma, u3, u4, stabit, sorted_eigs)

    #=A = fill(NaN, prod(dimvals), prod(dimvals))=#
    #=u1 = fill(NaN, dimvals[1], ranks[1])=#
    #=u2 = fill(NaN, dimvals[2], ranks[2])=#
    #=u3 = fill(NaN, p * dimvals[1], ranks[1])=#
    #=u4 = fill(NaN, p * dimvals[2], ranks[2])=#
    #=stabit = 0=#
    #==#
    #=u3_cont = Vector{Any}(undef, p)=#
    #=u4_cont = Vector{Any}(undef, p)=#
    #=while true=#
    #=    stabit += 1=#
    #=    randU1 = randn(dimvals[1], ranks[1])=#
    #=    u1 .= svd(randU1).U=#
    #=    randU2 = randn(dimvals[2], ranks[2])=#
    #=    u2 .= svd(randU2).U=#
    #=    u3 = randn(dimvals[1], ranks[1])=#
    #=    u3_scale = u3[1, 1]=#
    #=    u1 = copy(u1) * u3_scale=#
    #=    u3 = copy(u3) / u3_scale=#
    #=    u4 = randn(dimvals[2], ranks[2])=#
    #=    k21 = kron(u2, u1)[:, 1:prod(ranks)]=#
    #=    k43 = kron(u4, u3)[:, 1:prod(ranks)]=#
    #=    A .= k21 * k43'=#
    #=    if isstable(A, maxeigen)=#
    #=        break=#
    #=    end=#
    #=end=#
    #==#
    #=sorted_eigs = sort(abs.(eigvals(A)), rev=true)=#
    #=return (; A, u1, u2, u3, u4, stabit, sorted_eigs)=#
end

function simulate_rrmar_data(
    dimvals::AbstractVector,
    ranks::AbstractVector,
    obs::Int;
    A=nothing,
    snr=0.7,
    burnin=100,
    matrix_err=false,
)
    if isnothing(A)
        A = generate_rrmar_coef(dimvals, ranks)
    end

    dim1, dim2 = dimvals[1], dimvals[2]
    coef = A.A
    rho = maximum(abs.(eigvals(coef)))
    sigma_val = rho / snr
    diagerr = repeat([sigma_val], dim1 * dim2)
    sigma = diagm(diagerr)
    pre_sigma2, pre_sigma1 = nearest_kron(sigma, (dim2, dim2), (dim1, dim1))
    sigma1 = abs.(pre_sigma1) / norm(pre_sigma1)
    sigma2 = abs.(pre_sigma2) * norm(pre_sigma1)
    pre_data = zeros(dimvals[1] * dimvals[2], obs)
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
        for i = 2:obs
            vec_epsilon = rand(d)
            pre_data[:, i] .= coef * pre_data[:, i-1] + vec_epsilon
        end
        data = pre_data[:, (burnin+1):end]
        sigma = copy(diagerr)
        return (; data, coef, sigma)
    end
end
