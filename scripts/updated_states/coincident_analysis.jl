using DrWatson
@quickactivate :PseudoStructuralComovements
using CSV, Tables, SparseArrays, Distributions, XLSX
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

state_names = ["IA", "IL", "IN", "MI", "MN", "ND", "OH", "SD", "WI"]
rawdata = XLSX.readdata(datadir("./state_indexes/reguib_northcentral.xlsx"), "Sheet1!A2:S459")
vecdata = Float64.(rawdata[:, 2:end])'
coincident = vecdata[1:2:end, 278:end-2]'

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
# iowa, illinois, indiana, michigan, minnesota, north dakota, ohio, south dakota, wisconsin
rearranged_wages = ut_wages[:, [3,1,2,4,5,6,7,8,9]]
# monthly_wages = denton(rearranged_wages, ut_employment)
monthly_wages = quarterly_to_monthly(rearranged_wages)
wages = transform(monthly_wages)

catted_data = cat(employment', unemployment', hours', wages'; dims = 3)
tendata = permutedims(catted_data, (1,3,2))
# Corresponds to the perm WI, ND, OH, MN, MI, IA, SD, IL, IN
# perm_states = [9, 6, 7, 5, 4, 1, 8, 2, 3]
# Corresponds to the perm ND, IN, MI, MN, IL, SD, WI, OH, IA
perm_states = [6, 3, 4, 5, 2, 8, 9, 7, 1]
rearranged_tendata = tendata[perm_states, :, :]
n1, n2, obs = size(rearranged_tendata)
dimvals = [n1, n2]
matdata = vectorize(rearranged_tendata)

icest = rank_selection(matdata, dimvals; iters=1000, pmax=1)
icest2 = rrmar_ic(matdata, dimvals)
# aic selects (8,4) with one lag
# bic selects (5,4) with one lag
# ebic selects (2,1) with one lag
# this is 100 total parameters to estimate

n1, n2, obs = size(tendata)
dimvals = [n1, n2]
res = comovement_reg(matdata, dimvals, [2, 1]; iters=1000, p=1)
# save(datadir("updated_states/coincident_results.jld2"), Dict("res" => res))
# load the results
loaded_results = load(datadir("updated_states/coincident_results.jld2"))
res = loaded_results["res"]

# julia> res.delta_est
# 9×7 Matrix{Float64}:
#   1.0        0.0        0.0        0.0         0.0        0.0         0.0
#   0.0        1.0        0.0        0.0         0.0        0.0         0.0
#   0.0        0.0        1.0        0.0         0.0        0.0         0.0
#   0.0        0.0        0.0        1.0         0.0        0.0         0.0
#   0.0        0.0        0.0        0.0         1.0        0.0         0.0
#   0.0        0.0        0.0        0.0         0.0        1.0         0.0
#   0.0        0.0        0.0        0.0         0.0        0.0         1.0
#  -0.117196  -0.106524  -0.232348  -0.0724497  -0.403265   0.0258083  -0.0287248
#  -0.312692  -1.65608   -1.94897   -0.929204   -0.99046   -0.692681   -1.54464


# julia> res.gamma_est
# 4×3 Matrix{Float64}:
#   1.0        0.0      0.0
#   0.0        1.0      0.0
#   0.0        0.0      1.0
#  -0.632952  70.5833  -0.161184

# order is employment, unemployment, hours, wages
# employment co-moves with wages
# unemployment co-moves with wages
# hours do not co-move with wages

################################################################################
# Plots
################################################################################

# cis order is
# SD, MN, ND, IA, OH, IL, IN
# coincident order is 
# IA, IL, IN, MI, MN, ND, OH, SD, WI
# we omit MI and WI
factors = kron(res.u4_est, res.u3_est)' * matdata
# factors = kron(res.u4_est, I(4))' * matdata
u1 = nullspace(res.delta_est') * inv(nullspace(res.delta_est')[8:end, :])
# u2 = nullspace(res.gamma_est') * inv(nullspace(res.gamma_est')[4, 1])
# coincident_inds = kron(u2, u1) * factors
coincident_inds = u1 * factors
permed_cis = coincident_inds[invperm(perm_states), :]
# permed_cis = factors[invperm(perm_states), :]

save(datadir("updated_states/coincident_series.jld2"), Dict("cis" => permed_cis))

# Iowa
Plots.plot(demean_standardize(coincident[:, 1]), label = "Crone and Clayton-Matthews", title = "Iowa")
Plots.plot!(demean_standardize(permed_cis[1, :]), label = "Pseudo-Structural")
cor(coincident[:, 1], permed_cis[1, :])

# Illinois
Plots.plot(demean_standardize(coincident[:, 2]), label = "Crone and Clayton-Matthews", title = "Illinois")
Plots.plot!(demean_standardize(permed_cis[2, :]), label = "Pseudo-Structural")
cor(coincident[:, 2], permed_cis[2, :])

# Indiana
Plots.plot(demean_standardize(coincident[:, 3]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[3, :]), label = "Pseudo-Structural", title = "Indiana")
cor(coincident[:, 3], permed_cis[3, :])

# Michigan
Plots.plot(demean_standardize(coincident[:, 4]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[4, :]), label = "Pseudo-Structural", title = "Minnesota")
cor(coincident[:, 4], permed_cis[4, :])

# Minnesota
Plots.plot(demean_standardize(coincident[:, 5]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[5, :]), label = "Pseudo-Structural", title = "Minnesota")
cor(coincident[:, 5], permed_cis[5, :])

# North Dakota
Plots.plot(demean_standardize(coincident[:, 6]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[6, :]), label = "Pseudo-Structural", title = "North Dakota")
cor(coincident[:, 6], permed_cis[6, :])
# Plots.plot!(demean_standardize(employment[:, 5]), label = "Employment")
# ND experienced an oil boom between 2006, peaked in 2012, and crashed in 2014.
# can plot it as
# Plots.plot(ut_employment[:, 5])

# Ohio
Plots.plot(demean_standardize(coincident[:, 7]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[7, :]), label = "Pseudo-Structural", title = "Ohio")
cor(coincident[:, 7], permed_cis[7, :])

# South Dakota
Plots.plot(demean_standardize(coincident[:, 8]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[8, :]), label = "Pseudo-Structural", title = "South Dakota")
cor(coincident[:, 8], permed_cis[8, :])

# Wisconsin
Plots.plot(demean_standardize(coincident[:, 9]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[9, :]), label = "Pseudo-Structural", title = "Wisconsin")
cor(coincident[:, 9], permed_cis[9, :])

# Delta method using the asymptotic variance for delta and gamma
dg_var = inv(res.hess_est)[1:17, 1:17]

# I need delta star and gamma star
delta_star = res.delta_est[8:end, :]
gamma_star = res.gamma_est[end, :]'

# I want the asymptotic variance of kron(gamma_star, delta_star)'
m, n = size(reshape(gamma_star', :, 1))
p, q = size(delta_star')
commutation_mat = comm_mat(q, m)
large_comm = kron(kron(I(n), commutation_mat), I(p))
partial_gamma = kron(I(m * n), vec(delta_star'))
partial_delta = kron(vec(gamma_star'), I(p * q))
J = large_comm * hcat(partial_delta, partial_gamma)

# Delta result
delta_method_variance = J * dg_var * J'
delta_method_stderr = sqrt.(diag(delta_method_variance))

theta = vec(kron(gamma_star', delta_star'))
tstats = theta ./ delta_method_stderr

pvals = 2 .* (1 .- cdf.(Normal(), abs.(tstats)))

sig_5 = abs.(tstats) .> 1.372
hcat(theta, delta_method_stderr, tstats, pvals, sig_5)
# It is really only the second variable in the list that is significant
pvals_mat = reshape(pvals, p * m, n * q)

################################################################################

# Coincident and leading indicators
leading_fixed = leading[:, [2,3,6,7]]

coincident_and_leading = vcat(factors, leading_fixed')
n1, _, obs = size(rearranged_tendata)
dimvals = [n1, 2]
icest = rank_selection(coincident_and_leading, dimvals; iters=1000, pmax=1)

# julia> icest.ictable'
# 14×6 adjoint(::Matrix{Float64}) with eltype Float64:
#  -27497.0     -27342.9     -27434.8     1.0  1.0  1.0
#  -26151.4     -25959.6     -26074.0     2.0  1.0  1.0
#  -25457.9     -25235.3     -25368.1     3.0  1.0  1.0
#  -24686.7     -24440.1     -24587.2     4.0  1.0  1.0
#  -29175.5     -28911.8     -29069.1     5.0  1.0  1.0
#       2.0e12       2.0e12       2.0e12  6.0  1.0  1.0
#       2.0e12       2.0e12       2.0e12  7.0  1.0  1.0
#       2.0e12       2.0e12       2.0e12  1.0  2.0  1.0
#       2.0e12       2.0e12       2.0e12  2.0  2.0  1.0
#       2.0e12       2.0e12       2.0e12  3.0  2.0  1.0
#       2.0e12       2.0e12       2.0e12  4.0  2.0  1.0
#       2.0e12       2.0e12       2.0e12  5.0  2.0  1.0
#       2.0e12       2.0e12       2.0e12  6.0  2.0  1.0
#       2.0e12       2.0e12       2.0e12  7.0  2.0  1.0


# julia> icest.ictable'
# 36×6 adjoint(::Matrix{Float64}) with eltype Float64:
#  -52421.9  -52176.5  -51900.5  1.0  1.0  1.0
#  -52608.0  -52314.7  -51985.0  2.0  1.0  1.0
#  -52646.7  -52312.0  -51935.7  3.0  1.0  1.0
#  -52666.3  -52296.6  -51880.9  4.0  1.0  1.0
#  -52664.6  -52266.1  -51818.2  5.0  1.0  1.0
#  -52663.4  -52242.7  -51769.7  6.0  1.0  1.0
#  -52655.2  -52218.5  -51727.6  7.0  1.0  1.0
#  -52650.3  -52204.1  -51702.4  8.0  1.0  1.0
#  -52649.3  -52199.8  -51694.6  9.0  1.0  1.0
#  -52436.3  -52175.0  -51881.1  1.0  2.0  1.0
#  -52615.4  -52306.2  -51958.6  2.0  2.0  1.0
#  -52664.1  -52313.4  -51919.3  3.0  2.0  1.0
#  -52707.5  -52321.8  -51888.2  4.0  2.0  1.0
#  -52744.9  -52330.5  -51864.7  5.0  2.0  1.0
#  -52755.1  -52318.5  -51827.5  6.0  2.0  1.0
#  -52753.4  -52300.8  -51791.9  7.0  2.0  1.0
#  -52755.1  -52293.0  -51773.4  8.0  2.0  1.0
#  -52755.3  -52289.9  -51766.8  9.0  2.0  1.0
#  -52440.3  -52169.4  -51864.8  1.0  3.0  1.0
#  -52623.3  -52304.5  -51946.2  2.0  3.0  1.0
#  -52685.2  -52325.0  -51920.1  3.0  3.0  1.0
#  -52736.5  -52341.2  -51896.9  4.0  3.0  1.0
#  -52780.7  -52356.8  -51880.2  5.0  3.0  1.0
#  -52795.0  -52348.7  -51847.0  6.0  3.0  1.0
#  -52805.7  -52343.5  -51823.9  7.0  3.0  1.0
#  -52807.8  -52336.1  -51805.8  8.0  3.0  1.0
#  -52807.9  -52333.0  -51799.0  9.0  3.0  1.0
#  -52439.8  -52165.6  -51857.5  1.0  4.0  1.0
#  -52621.7  -52299.8  -51937.8  2.0  4.0  1.0
#  -52700.5  -52337.2  -51928.6  3.0  4.0  1.0
#  -52748.6  -52350.1  -51902.2  4.0  4.0  1.0
#  -52794.3  -52367.2  -51887.0  5.0  4.0  1.0
#  -52805.3  -52355.9  -51850.6  6.0  4.0  1.0
#  -52812.9  -52347.5  -51824.3  7.0  4.0  1.0
#  -52813.9  -52339.0  -51805.0  8.0  4.0  1.0
#  -52813.4  -52335.2  -51797.7  9.0  4.0  1.0


# julia> icest2.ictable'
# 36×6 adjoint(::Matrix{Float64}) with eltype Float64:
#  -52415.0  -52169.6  -52405.6  1.0  1.0  1.0
#  -52601.1  -52307.9  -52523.1  2.0  1.0  1.0
#  -52639.8  -52305.1  -52502.3  3.0  1.0  1.0
#  -52659.4  -52289.7  -52471.6  4.0  1.0  1.0
#  -52657.7  -52259.3  -52428.7  5.0  1.0  1.0
#  -52656.5  -52235.8  -52395.5  6.0  1.0  1.0
#  -52648.2  -52211.5  -52364.3  7.0  1.0  1.0
#  -52643.2  -52197.0  -52345.7  8.0  1.0  1.0
#  -52642.2  -52192.8  -52340.1  9.0  1.0  1.0
#  -52426.8  -52165.4  -52390.4  1.0  2.0  1.0
#  -52605.5  -52296.3  -52500.5  2.0  2.0  1.0
#  -52655.5  -52304.8  -52491.0  3.0  2.0  1.0
#  -52700.6  -52315.0  -52485.9  4.0  2.0  1.0
#  -52738.1  -52323.7  -52482.2  5.0  2.0  1.0
#  -52748.4  -52311.7  -52460.4  6.0  2.0  1.0
#  -52745.6  -52293.0  -52434.8  7.0  2.0  1.0
#  -52748.3  -52286.1  -52423.8  8.0  2.0  1.0
#  -52748.5  -52283.1  -52419.4  9.0  2.0  1.0
#  -52430.3  -52159.4  -52377.8  1.0  3.0  1.0
#  -52613.6  -52294.8  -52492.5  2.0  3.0  1.0
#  -52673.1  -52312.9  -52492.5  3.0  3.0  1.0
#  -52723.3  -52328.0  -52492.4  4.0  3.0  1.0
#  -52767.5  -52343.6  -52495.4  5.0  3.0  1.0
#  -52782.4  -52336.2  -52478.4  6.0  3.0  1.0
#  -52792.7  -52330.5  -52465.7  7.0  3.0  1.0
#  -52795.4  -52323.6  -52454.7  8.0  3.0  1.0
#  -52795.8  -52320.9  -52450.6  9.0  3.0  1.0
#  -52430.1  -52155.9  -52372.1  1.0  4.0  1.0
#  -52612.3  -52290.4  -52485.8  2.0  4.0  1.0
#  -52689.9  -52326.5  -52503.9  3.0  4.0  1.0
#  -52737.3  -52338.9  -52501.0  4.0  4.0  1.0
#  -52783.3  -52356.2  -52505.9  5.0  4.0  1.0
#  -52794.3  -52344.9  -52484.9  6.0  4.0  1.0
#  -52800.9  -52335.5  -52468.6  7.0  4.0  1.0
#  -52801.9  -52327.0  -52455.8  8.0  4.0  1.0
#  -52801.3  -52323.2  -52450.7  9.0  4.0  1.0





