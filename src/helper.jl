
function vec_to_ll(v, n)
    @assert length(v) == n * (n + 1) / 2 "Vector length is wrong!"
    mat = zeros(eltype(v), n, n)
    idx = 1
    for j = 1:n
        num_elements = n - j + 1
        mat[j:end, j] = v[idx:(idx+num_elements-1)]
        idx += num_elements
    end
    return mat
end

function vech(m)
    n = size(m, 1)
    return vcat([m[j:end, j] for j in 1:n]...)
end

function rotate_u!(U)
    col_num = size(U, 2)

    U .= U * inv(U[1:col_num, 1:col_num])
    return U
end

function make_companion(omega, pi_mat; ll=nothing)
    k = size(pi_mat, 1)
    p = size(pi_mat, 2) ÷ size(pi_mat, 1)
    omega_tilde = diagm(ones(k * p))
    omega_tilde[1:k, 1:k] .= omega
    pi_tilde = diagm(-k => ones(k * p - k))
    pi_tilde[1:k, :] .= pi_mat
    ll_tilde = nothing
    if !isnothing(ll)
        ll_tilde = diagm(ones(k * p))
        ll_tilde[1:k, 1:k] .= ll
    end
    return (; omega_tilde, pi_tilde, ll_tilde)
end

function companion_data(data, p)
    p_adj = p - 1
    k, obs = size(data)
    obs_eff = obs - p_adj
    Z = Matrix{Float64}(undef, k * (p_adj + 1), obs_eff)

    for i in 0:p_adj
        Z[(i*k+1):(i+1)*k, :] = data[:, (p_adj-i+1):end-i]
    end
    return Z
end

function insertk!(v::AbstractVector, k::Int; val=1)
    n = length(v)
    starts = collect(1:k:n)
    for s in reverse(starts)
        insert!(v, s, val)
    end
    return v
end

function removek!(v::AbstractVector, k::Int; val=1)
    step = k + 1
    idxs = collect(1:step:length(v))
    for i in reverse(idxs)
        if v[i] == val
            deleteat!(v, i)
        end
    end
    return v
end

function vecb(M::AbstractMatrix, b::Integer)
    n, p = size(M)
    @assert n % b == 0 "number of rows ($(n)) must be a multiple of block‐size b=$(b)"
    nb = n ÷ b
    out = similar(vec(M))
    pos = 1
    for i in 0:(nb-1)
        rows = (i*b+1):(i*b+b)
        for j in 1:p
            @views out[pos:pos+b-1] .= M[rows, j]
            pos += b
        end
    end
    return out
end

function cov_to_ll(mat)
    N = size(mat, 1)
    L = cholesky(mat).L
    return [L[i, j] for j = 1:N for i = j:N]
end

function ll_to_cov(vecL::AbstractVector, n)
    L = vec_to_ll(vecL, n)
    return L * L'
end

function check_conf(v, lower, upper)
    return Int.((v .>= lower) .& (v .<= upper))
end

function var_coef(data)
    cen_data = data .- mean(data, dims=2)
    resp = cen_data[:, 2:end]
    pred = cen_data[:, 1:(end-1)]

    coef = resp * pred' * inv(pred * pred')
    return coef
end

function ols_coef(resp, pred)
    cen_resp = resp .- mean(resp, dims=2)
    cen_pred = pred .- mean(pred, dims=2)

    coef = cen_resp * cen_pred' * inv(cen_pred * cen_pred')
    return coef
end

function nearest_posdef(A; minimumeigenvalue=1e-6)
    Evalues, Evector = eigen(A)
    newA = Evector * Diagonal(max.(Evalues, minimumeigenvalue)) * Evector'

    return Symmetric(newA)
end

function nearest_kron(A, SizeB, SizeC)
    m, n = size(A)
    m1, n1 = SizeB
    m2, n2 = SizeC

    @assert m == m1 * m2 "Size mismatch: m ≠ m1 * m2"
    @assert n == n1 * n2 "Size mismatch: n ≠ n1 * n2"

    ten_A = reshape(A, (m2, m1, n2, n1))
    R = reshape(permutedims(ten_A, (2, 4, 1, 3)), m1 * n1, m2 * n2)

    F = svd(R)
    B = reshape(F.U[:, 1] * sqrt(F.S[1]), m1, n1)
    C = reshape(F.V[:, 1] * sqrt(F.S[1]), m2, n2)
    Aest = kron(B, C)
    norm_a = norm(A - Aest)
    return B, C, norm_a
end
