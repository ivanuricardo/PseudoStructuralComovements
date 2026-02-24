
function load_series(kind::AbstractString; col=5)
    dir = datadir("updated_states/")
    files = filter(name -> endswith(name, "-$(kind).csv"), readdir(dir))
    n = length(files)

    col_vectors = Vector{Vector{Float64}}(undef, n)

    for (i, f) in enumerate(files)
        tbl = Tables.columntable(CSV.File(joinpath(dir, f)))
        v = tbl[col]
        col_vectors[i] = Float64.(v)
    end

    maxlen = maximum(length.(col_vectors))
    X = fill(NaN, maxlen, n)

    for (j, v) in enumerate(col_vectors)
        X[1:length(v), j] .= v
    end

    return X
end

function transform(data::AbstractMatrix; type::String = "logdiff")
    N, P = size(data)
    N > 1 || error("Data must have at least 2 rows")

    if type == "diff"
        return data[2:end, :] .- data[1:end-1, :]

    elseif type == "logdiff"
        any(data .<= 0) && error("Log differences require strictly positive data")
        return log.(data[2:end, :]) .- log.(data[1:end-1, :])

    else
        error("Unknown type: $type. Use \"diff\" or \"logdiff\".")
    end
end

function monthly_to_quarterly(X::AbstractMatrix; agg_type::String="mean")
    n, m = size(X)
    if n % 3 != 0
        error("Number of rows ($n) is not divisible by 3. Cannot form full quarters.")
    end
    n_quarters = div(n, 3)
    Q = zeros(n_quarters, m)
    @inbounds for i in 1:n_quarters
        rows = (3*(i-1)+1):(3*i)
        if agg_type == "mean"
            Q[i, :] = mean(X[rows, :], dims=1)
        elseif agg_type == "sum"
            Q[i, :] = sum(X[rows, :], dims=1)
        else
            error("Unknown aggregation type: $agg_type. Use 'mean' or 'sum'.")
        end
    end
    return Q
end

function vectorize(data::AbstractArray{T}) where T
    n1, n2, obs = size(data)
    return reshape(data, n1*n2, obs)
end

function demean_standardize(X; dims=1, corrected=true)
    μ = mean(X; dims=dims)
    σ = std(X; dims=dims, corrected=corrected)

    return (X .- μ) ./ σ
end

function denton(Q::AbstractVector, x::AbstractVector; s::Int=3)
    T = length(Q)
    months = length(x)
    @assert months == s * T "length(x) must equal s * length(Q)"

    # Ensure positivity for proportional Denton
    xv = Float64.(x)
    @assert all(xv .>= 0.0) "indicator contains nonpositive values"

    # Aggregation matrix J (T x months)
    J = Float64.(kron(I(T), ones(s)'))  # size (n x months)

    # First-difference matrix D ((months-1) x months)
    D = zeros(Float64, months-1, months)
    for i in 1:(months-1)
        D[i,i] = -1.0
        D[i,i+1] =  1.0
    end
    D = sparse(D)

    # Proportional weighting matrix
    W = spdiagm(0 => 1.0 ./ xv)

    # sandwiched matrix
    M = (D * W)' * (D * W)

    # KKT system
    KKT = [M  J';
           J  spzeros(T,T)]
    rhs = vcat(zeros(months), Q)

    sol = KKT \ rhs
    return sol[1:months]
end

function denton(Q::AbstractMatrix, x::AbstractMatrix; s::Int=3)
    T, n_series = size(Q, 1), size(Q, 2)
    months = size(x, 1)
    @assert months == s * T "length(x) must equal s * size(Q,1)"
    @assert size(x,2) == n_series "x and Q must have the same number of columns"

    Y = zeros(months, n_series)
    for i in 1:n_series
        Y[:, i] = denton(Q[:, i], x[:, i]; s=s)
    end
    return Y
end

function quarterly_to_monthly(A::AbstractMatrix)
    T, p = size(A)
    Q = A ./ 3
    T3 = 3 * T
    R = Array{eltype(Q)}(undef, T3, p)
    @inbounds for i in 1:T
        rows = (3*(i-1)+1):(3*i)
        R[rows, :] .= repeat(Q[i, :]', 3, 1)
    end
    return R
end

function comm_mat(m::Int, n::Int)
    A = reshape(1:m*n, m, n)
    v = vec(A')

    P = Matrix{Float64}(I, m*n, m*n)
    return P[v, :]
end

function monthly_aggregation(kind::AbstractString)
    dir = datadir("updated_states/")
    files = filter(name -> endswith(name, "-$(kind).csv"), readdir(dir))
    n = length(files)
    col_vectors = Vector{Vector{Float64}}(undef, n)

    for (i, f) in enumerate(files)
        df = CSV.read(joinpath(dir, f), DataFrame)
        col_sym = names(df)[2]
        df.month = Date.(year.(df.observation_date), month.(df.observation_date), 1)
        monthly = combine(groupby(df, :month), col_sym => sum => :monthly_claims)
        col_vectors[i] = monthly.monthly_claims
    end

    maxlen = maximum(length.(col_vectors))
    X = fill(NaN, maxlen, n)

    for (j, v) in enumerate(col_vectors)
        X[1:length(v), j] .= v
    end

    return X
end

