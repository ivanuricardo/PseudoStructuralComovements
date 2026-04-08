using LinearAlgebra
using Statistics

struct RRVARFit
    p::Int
    r::Int
    intercept::Bool
    α::Vector{Float64}
    B::Matrix{Float64}
    Σ::Matrix{Float64}
    loglik::Float64
    aic::Float64
    bic::Float64
    ebic::Float64
end

function rrvar_num_params(n, p, r; intercept=true)
    k = n * p
    nB = r * (k + n - r)
    nα = intercept ? n : 0
    nΣ = n * (n + 1) ÷ 2
    return nB + nα + nΣ
end

"Build response matrix Ydep and lagged regressor matrix X for a VAR(p)."
function lagmat(Y::AbstractMatrix, p::Int)
    T, n = size(Y)
    @assert p ≥ 1 "p must be at least 1"
    @assert T > p "Need T > p"

    Ydep = Matrix{Float64}(Y[p+1:end, :])

    blocks = [Matrix{Float64}(Y[p+1-j:end-j, :]) for j in 1:p]
    X = hcat(blocks...)

    return Ydep, X
end

"Symmetric positive definite inverse square root."
function invsqrt_spd(S::AbstractMatrix)
    F = eigen(Symmetric(Matrix(S)))
    @assert all(F.values .> 0) "Matrix must be positive definite"
    return F.vectors * Diagonal(1.0 ./ sqrt.(F.values)) * F.vectors'
end

"Symmetric positive definite square root."
function sqrt_spd(S::AbstractMatrix)
    F = eigen(Symmetric(Matrix(S)))
    @assert all(F.values .> 0) "Matrix must be positive definite"
    return F.vectors * Diagonal(sqrt.(F.values)) * F.vectors'
end

"""
Reduced-rank regression estimator:
min ||Y - X*B||_F^2 subject to rank(B) ≤ r

Returns B (k×m).
"""
function reduced_rank_regression(X::AbstractMatrix, Y::AbstractMatrix, r::Int)
    T = size(X, 1)

    Sxx = (X' * X) / T
    Syy = (Y' * Y) / T
    Syx = (Y' * X) / T

    Wx = invsqrt_spd(Sxx)
    Wy = invsqrt_spd(Syy)

    M = Wy * Syx * Wx
    F = svd(M)

    r = min(r, length(F.S))
    U = F.U[:, 1:r]
    V = F.Vt'[:, 1:r]
    D = Diagonal(F.S[1:r])

    B = Wx * V * D * U' * sqrt_spd(Syy)
    return B
end

"Gaussian log-likelihood given residual matrix E and covariance Σ."
function gaussian_loglik(E::AbstractMatrix, Σ::AbstractMatrix)
    T, n = size(E)
    C = cholesky(Symmetric(Matrix(Σ)))
    logdetΣ = 2sum(log, diag(C.U))

    quad = 0.0
    for t in 1:T
        e = vec(E[t, :])
        quad += dot(e, C \ e)
    end

    return -0.5 * (T * n * log(2π) + T * logdetΣ + quad)
end

"""
Fit a reduced-rank VAR(p).

If intercept=true, the intercept is estimated unrestricted after centering.
The rank restriction is applied only to the lag coefficient matrix.
"""
function fit_rrvar(Y::AbstractMatrix, p::Int, r::Int; intercept::Bool=true)
    Ydep, X = lagmat(Y, p)
    n = size(Ydep, 2)

    if intercept
        μY = vec(mean(Ydep, dims=1))
        μX = vec(mean(X, dims=1))

        Yc = Ydep .- μY'
        Xc = X .- μX'

        B = reduced_rank_regression(Xc, Yc, r)
        α = vec(μY' .- μX' * B)

        Ŷ = X * B .+ α'
    else
        B = reduced_rank_regression(X, Ydep, r)
        α = zeros(n)

        Ŷ = X * B
    end

    E = Ydep .- Ŷ
    Σ = (E' * E) / size(E, 1)
    ll = gaussian_loglik(E, Σ)

    T_eff = size(E, 1)
    n_params = rrvar_num_params(n, p, r; intercept=intercept)

    aic = -2 * ll + 2 * n_params
    bic = -2 * ll + n_params * log(T_eff)
    ebic = -2 * ll + n_params * log(T_eff) + log(n)

    return RRVARFit(p, r, intercept, α, B, Σ, ll, aic, bic, ebic)

end

"Log-likelihood for a fitted model on new data Y."
function rrvar_loglik(Y::AbstractMatrix, fit::RRVARFit)
    Ydep, X = lagmat(Y, fit.p)
    Ŷ = fit.intercept ? (X * fit.B .+ fit.α') : (X * fit.B)
    E = Ydep .- Ŷ
    return gaussian_loglik(E, fit.Σ)
end




Y = randn(300, 3)          # T × n data matrix, rows = time
fit = fit_rrvar(Y, 2, 1)   # VAR(2), rank 1

fit.B
fit.α
fit.Σ
fit.loglik
