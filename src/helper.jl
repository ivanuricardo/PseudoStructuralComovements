
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

function cov_to_ll(mat)
    N = size(mat, 1)
    L = cholesky(mat).L
    return [L[i, j] for j = 1:N for i = j:N]
end

function ll_to_cov(vecL::Vector, n::Int)
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
