using DrWatson
@quickactivate :PseudoStructuralComovements
using CSV, Tables, SparseArrays, XLSX
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

rawdata = XLSX.readdata(datadir("./state_indexes/reguib_northcentral.xlsx"), "Sheet1!A2:S459")
vecdata = Float64.(rawdata[:, 2:end])'
coincident = vecdata[1:2:end, 232:end]'

# Before transformation
ut_employment = load_series("employment")
ut_unemployment = load_series("unemployment")
ut_hours = load_series("hours")

employment = transform(ut_employment)
unemployment = transform(ut_unemployment; type = "diff")
hours = transform(ut_hours)

file = XLSX.readxlsx(datadir("updated_states/wages.xlsx"))
sheet = file["Table"]
ut_wages = Float64.(XLSX.getdata(sheet, "C7:BJ15"))'
# Rearrange to be alphabetical and align with other series
# ia, il, in, mi, mn, nd, oh, sd, wi
rearranged_wages = ut_wages[:, [3,1,2,4,5,6,7,8,9]]
# monthly_wages = denton(rearranged_wages, ut_employment)
monthly_wages = quarterly_to_monthly(rearranged_wages)
wages = transform(monthly_wages)

ser = 1:9
Plots.plot(demean_standardize(ut_employment[:, ser]))
Plots.plot(demean_standardize(unemployment[:, ser]))
Plots.plot!(demean_standardize(wages[:, ser]))
Plots.plot(ut_hours[:, 3])
Plots.plot(coincident[:, ser])

catted_data = cat(wages', unemployment', hours', employment'; dims = 3)
tendata = permutedims(catted_data, (1,3,2))
# Corresponds to the perm WI, MI, SD, MN, ND, IA, OH, IL, IN
pp = [9, 4, 8, 5, 6, 1, 7, 2, 3]
# pp = [2,3,6]
# or perm can be IA, IL, MN, OH, SD, ND, IL
# pp = [1,2,4,6,7,3,5]
rearranged_tendata = tendata[pp, :, :]
matdata = vectorize(rearranged_tendata)

n1, n2, obs = size(rearranged_tendata)
dimvals = [n1, n2]
icest = rank_selection(matdata, dimvals; iters=1000, pmax=1)
icest = rrmar_ic(matdata, dimvals)
# aic selects (8,4) with one lag
# bic selects (5,4) with one lag
# ebic selects (2,1) with one lag
# this is 100 total parameters to estimate

