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

ut_claims = seasonal_adjust(seas_claims)
claims = transform(ut_claims)

Plots.plot(demean_standardize(claims))
Plots.plot(claims[:, 1])


