using DrWatson
@quickactivate :PseudoStructuralComovements
using CSV, Tables
using XLSX
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

# Without using Michigan or Wisconsin
# Not sure if I should still keep this series
rawdata = XLSX.readdata(datadir("./state_indexes/reguib_northcentral.xlsx"), "Sheet1!A2:S459")
vecdata = Float64.(rawdata[:, 2:end])'
untransformed_coincident = vecdata[1:2:end, 232:end-2]'
quarterly_coincident = monthly_to_quarterly(untransformed_coincident)
coincident = quarterly_coincident[:, [1,2,3,5,6,7,8]]

file = XLSX.readxlsx(datadir("updated_states/wages.xlsx"))
sheet = file["Table"]
ut_wages = Float64.(XLSX.getdata(sheet, "O7:CL13"))'
# Rearrange to be alphabetical and align with other series
# iowa, illinois, indiana, minnesotta, north dakota, ohio, south dakota
# maybe I could try without north dakota
rearranged_wages = ut_wages[:, [3,1,2,4,5,6,7]]
wages = transform(rearranged_wages[:, [1,2,3,4,6,7]])

start_date = 133  # Jan 2001
covid_start = 12  # Dec 2019

# Before transformation
ut_employment = load_series("employment")[start_date:(end-covid_start), :]
ut_unemployment = load_series("unemployment")[start_date:(end-covid_start), :]
ut_hours = load_series("hours")[1:(end-covid_start), :]  # already starts in Jan

# Aggregate to quarterly
quarterly_employment = monthly_to_quarterly(ut_employment)
quarterly_unemployment = monthly_to_quarterly(ut_unemployment)
quarterly_hours = monthly_to_quarterly(ut_hours)

employment = transform(quarterly_employment)[:, [1,2,3,4,6,7]]
unemployment = transform(quarterly_unemployment; type = "diff")[:, [1,2,3,4,6,7]]
hours = transform(quarterly_hours)[:, [1,2,3,4,6,7]]

catted_data = cat(employment', unemployment', hours', wages'; dims = 3)
tendata = permutedims(catted_data, (1,3,2))
matdata = vectorize(tendata)

ser = 1
Plots.plot(demean_standardize(employment[:, ser]))
Plots.plot!(-demean_standardize(unemployment[:, ser]))
Plots.plot!(demean_standardize(hours[:, ser]))
Plots.plot!(demean_standardize(wages[:, ser]))
Plots.plot!(demean_standardize(coincident[:, ser]))

n1, n2, obs = size(tendata)
dimvals = [n1, n2]
icest = rank_selection(matdata, dimvals; iters=1000, pmax=1)
# aic selects (6,4) with one lag
# bic selects (6,4) with one lag (5,4 without ND)
# hqc selects (6,4) with one lag
# this is 100 total parameters to estimate

ser = 5
Plots.plot(hcat(employment[:, ser], wages[:, ser]))

res = comovement_reg(matdata, dimvals, [2, 1]; iters=1000, p=1)
# Nothing is significant except the relation between Iowa and Ohio
# I could e.g., restrict the second dimension to be 1. Then I get two significant coefs
# however, if I restrict the second dimension to be 1, then the first dimension
# should be 3 or 4 (3 for HQC 4 for AIC, 2 for BIC)
#
#  employment unemployment hours wages
# julia> res.gamma_est  # for rank (3,1)
# 4×3 Matrix{Float64}:
#   1.0        0.0     0.0
#   0.0        1.0     0.0
#   0.0        0.0     1.0
#  -0.683243  75.431  -0.207048
# julia> res.gamma_est  # for rank (2,1)
# 4×3 Matrix{Float64}:
#   1.0        0.0      0.0
#   0.0        1.0      0.0
#   0.0        0.0      1.0
#  -0.724399  66.9276  -0.125401
#
# julia> res.gamma_est  # for rank (7,1)
# 4×3 Matrix{Float64}:
#   1.0        0.0      0.0
#   0.0        1.0      0.0
#   0.0        0.0      1.0
#  -0.845962  99.0915  -0.214233
# The estimates are quite similar, as well as std errors.
#
# julia> icest.ictable'
# 28×6 adjoint(::Matrix{Float64}) with eltype Float64:
#   AIC       BIC       HQC       ranks    lag
#  -15094.5  -14964.7  -15042.7  1.0  1.0  1.0
#  -15146.0  -14990.7  -15084.0  2.0  1.0  1.0
#  -15160.9  -14984.7  -15090.5  3.0  1.0  1.0
#  -15164.6  -14972.2  -15087.8  4.0  1.0  1.0
#  -15164.5  -14960.5  -15083.1  5.0  1.0  1.0
#  -15164.3  -14953.4  -15080.1  6.0  1.0  1.0
#  -15162.7  -14949.5  -15077.6  7.0  1.0  1.0
#  -15106.1  -14964.7  -15049.7  1.0  2.0  1.0
#  -15175.5  -15008.6  -15108.9  2.0  2.0  1.0
#  -15199.6  -15011.9  -15124.7  3.0  2.0  1.0
#  -15227.4  -15023.4  -15145.9  4.0  2.0  1.0
#  -15349.5  -15134.0  -15263.5  5.0  2.0  1.0
#  -15348.4  -15126.0  -15259.6  6.0  2.0  1.0
#  -15346.6  -15121.8  -15256.9  7.0  2.0  1.0
#  -15106.5  -14958.2  -15047.3  1.0  3.0  1.0
#  -15247.8  -15074.0  -15178.4  2.0  3.0  1.0
#  -15327.3  -15132.6  -15249.6  3.0  3.0  1.0
#  -15359.4  -15148.5  -15275.2  4.0  3.0  1.0
#  -15391.7  -15169.3  -15302.9  5.0  3.0  1.0
#  -15405.3  -15175.8  -15313.7  6.0  3.0  1.0
#  -15405.3  -15173.5  -15312.7  7.0  3.0  1.0
#  -15105.6  -14954.9  -15045.4  1.0  4.0  1.0
#  -15249.8  -15073.7  -15179.5  2.0  4.0  1.0
#  -15331.5  -15134.5  -15252.8  3.0  4.0  1.0
#  -15364.3  -15151.0  -15279.1  4.0  4.0  1.0
#  -15398.0  -15173.2  -15308.3  5.0  4.0  1.0
#  -15415.5  -15183.8  -15323.0  6.0  4.0  1.0
#  -15415.1  -15181.0  -15321.6  7.0  4.0  1.0




