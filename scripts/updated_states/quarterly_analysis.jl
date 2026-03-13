using DrWatson
@quickactivate :PseudoStructuralComovements
using CSV, Tables
using XLSX
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

rawdata = XLSX.readdata(datadir("./state_indexes/reguib_northcentral.xlsx"), "Sheet1!A2:S459")
vecdata = Float64.(rawdata[:, 2:end])'
untransformed_coincident = vecdata[1:2:end, 277:end-2]'
quarterly_coincident = monthly_to_quarterly(untransformed_coincident)
coincident = quarterly_coincident[:, [1,2,3,5,6,7,8]]

file = XLSX.readxlsx(datadir("updated_states/wages.xlsx"))
sheet = file["Table"]
ut_wages = Float64.(XLSX.getdata(sheet, "C7:BJ15"))'
# Rearrange to be alphabetical and align with other series
# iowa, illinois, indiana, michigan, minnesota, north dakota, ohio, south dakota
rearranged_wages = ut_wages[:, [3,1,2,4,5,6,7,8,9]]
wages = transform(rearranged_wages)

# Before transformation
ut_employment = load_series("employment")
ut_unemployment = load_series("unemployment")
ut_hours = load_series("hours")

# Aggregate to quarterly
quarterly_employment = monthly_to_quarterly(ut_employment)
quarterly_unemployment = monthly_to_quarterly(ut_unemployment)
quarterly_hours = monthly_to_quarterly(ut_hours)

employment = transform(quarterly_employment)
unemployment = transform(quarterly_unemployment; type = "diff")
hours = transform(quarterly_hours)

catted_data = cat(employment', unemployment', hours', wages'; dims = 3)
tendata = permutedims(catted_data, (1,3,2))
matdata = vectorize(tendata)

ser = 4
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

# Just for  9 states, quarterly
# julia> icest.ictable'
# 36×6 adjoint(::Matrix{Float64}) with eltype Float64:
#  -15457.8  -15297.8  -15475.0  1.0  1.0  1.0
#  -15505.4  -15314.2  -15470.6  2.0  1.0  1.0
#  -15531.3  -15313.2  -15451.6  3.0  1.0  1.0
#  -15538.8  -15297.8  -15420.9  4.0  1.0  1.0
#  -15542.6  -15282.9  -15393.5  5.0  1.0  1.0
#  -15539.0  -15264.8  -15365.7  6.0  1.0  1.0
#  -15530.2  -15245.6  -15339.6  7.0  1.0  1.0
#  -15532.1  -15241.2  -15331.1  8.0  1.0  1.0
#  -15530.7  -15237.8  -15326.2  9.0  1.0  1.0
#  -15469.0  -15298.7  -15464.9  1.0  2.0  1.0
#  -15580.7  -15379.2  -15524.6  2.0  2.0  1.0
#  -15642.7  -15414.1  -15541.5  3.0  2.0  1.0
#  -15664.0  -15412.6  -15524.7  4.0  2.0  1.0
#  -15594.7  -15324.6  -15424.2  5.0  2.0  1.0
#  -15616.9  -15332.2  -15422.2  6.0  2.0  1.0
#  -15625.8  -15330.7  -15413.7  7.0  2.0  1.0
#  -15628.2  -15326.9  -15405.8  8.0  2.0  1.0
#  -15680.5  -15377.1  -15454.6  9.0  2.0  1.0
#  -15476.0  -15299.4  -15459.0  1.0  3.0  1.0
#  -15592.0  -15384.3  -15523.1  2.0  3.0  1.0
#  -15665.5  -15430.7  -15551.5  3.0  3.0  1.0
#  -15701.8  -15444.1  -15549.7  4.0  3.0  1.0
#  -15728.8  -15452.5  -15545.5  5.0  3.0  1.0
#  -15747.8  -15457.0  -15540.3  6.0  3.0  1.0
#  -15760.8  -15459.6  -15536.0  7.0  3.0  1.0
#  -15766.2  -15458.8  -15531.0  8.0  3.0  1.0
#  -15764.7  -15455.2  -15526.0  9.0  3.0  1.0
#  -15475.2  -15296.6  -15453.9  1.0  4.0  1.0
#  -15592.3  -15382.4  -15519.0  2.0  4.0  1.0
#  -15664.0  -15427.2  -15545.8  3.0  4.0  1.0
#  -15703.2  -15443.5  -15546.8  4.0  4.0  1.0
#  -15728.6  -15450.2  -15541.1  5.0  4.0  1.0
#  -15748.6  -15455.7  -15536.8  6.0  4.0  1.0
#  -15764.4  -15461.1  -15535.3  7.0  4.0  1.0
#  -15771.0  -15461.5  -15531.5  8.0  4.0  1.0
#  -15769.6  -15458.0  -15526.7  9.0  4.0  1.0
# BIC selects 8,4, ebic selects 3,3



