using DrWatson
@quickactivate :PseudoStructuralComovements
using Tables, SparseArrays
using CSV, DataFrames, Dates, RCall
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

seas_claims = monthly_aggregation("claims")

function seasonal_adjust(matrix::Matrix{Float64})
    @rput matrix
    R"""
    library(seasonal)
    start_year = 2001
    start_month = 1
    ts_list <- lapply(1:ncol(matrix), function(i) {
        ts(matrix[, i], start = c(start_year, start_month), frequency = 12)
    })
    sa_list <- lapply(ts_list, seas)
    sa_values <- sapply(sa_list, function(x) final(x))
    """
    @rget sa_values
    sa_values
end

# log differences
ut_claims = seasonal_adjust(seas_claims)
claims = transform(ut_claims)

# log differences
loaded_results = load(datadir("updated_states/permits.jld2"))
seas_permits = loaded_results["res"]
ut_permits = seasonal_adjust(seas_permits)
permits = transform(ut_permits)

# no transformation
yield_dir = datadir("updated_states/yield_minus_ffr.csv")
ut_yield_minus_ffr = CSV.read(yield_dir, DataFrame; header = 1)
yield_minus_ffr = repeat(ut_yield_minus_ffr.T10YFFM, 1, 7)

loaded_cis = load(datadir("updated_states/coincident_series.jld2"))
coincident_inds = loaded_cis["cis"]

li_data = vcat(coincident_inds, claims', permits', yield_minus_ffr[2:end, 1]')
save(datadir("updated_states/leading_data.jld2"), Dict("li_data" => li_data))
