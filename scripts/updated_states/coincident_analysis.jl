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

# catted_data = cat(wages', hours', unemployment', employment'; dims = 3)
catted_data = cat(employment', unemployment', hours', wages'; dims = 3)
tendata = permutedims(catted_data, (1,3,2))
# Corresponds to the perm WI, ND, OH, MN, MI, IA, SD, IL, IN
# perm_states = [9, 6, 7, 5, 4, 1, 8, 2, 3]
# Corresponds to the perm IA, IN, MI, MN, OH, SD, WI, IL, ND
perm_states = [1, 3, 4, 5, 7, 8, 9, 2, 6]
rearranged_tendata = tendata[perm_states, :, :]
matdata = vectorize(rearranged_tendata)

n1, n2, obs = size(tendata)
dimvals = [n1, n2]
res = comovement_reg(matdata, dimvals, [2, 1]; iters=1000, p=1)
# load the results
loaded_results = load(datadir("updated_states/coincident_results.jld2"))
other_res = loaded_results["res"]

# julia> res.delta_est
# 9×7 Matrix{Float64}:
#   1.0         0.0         0.0         0.0          0.0        0.0         0.0
#   0.0         1.0         0.0         0.0          0.0        0.0         0.0
#   0.0         0.0         1.0         0.0          0.0        0.0         0.0
#   0.0         0.0         0.0         1.0          0.0        0.0         0.0
#   0.0         0.0         0.0         0.0          1.0        0.0         0.0
#   0.0         0.0         0.0         0.0          0.0        1.0         0.0
#   0.0         0.0         0.0         0.0          0.0        0.0         1.0
#   0.0823831  -0.0923376   0.0757538  -0.00634231  -0.110979   0.0709135  -0.338586
#  -0.942411   -1.18696    -0.423185   -0.566303    -0.172774  -0.61059    -0.601133


# julia> res.gamma_est
# 4×3 Matrix{Float64}:
#   1.0        0.0     0.0
#   0.0        1.0     0.0
#   0.0        0.0     1.0
#  -1.68575  105.162  -0.348001

# order is wages, unemployment, hours, employment
# wages co-moves with employment
# unemployment co-moves with employment
# hours do not co-move with employment

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
Plots.plot!(-demean_standardize(permed_cis[1, :]), label = "Pseudo-Structural")
cor(coincident[:, 1], -permed_cis[1, :])

# Illinois
Plots.plot(demean_standardize(coincident[:, 2]), label = "Crone and Clayton-Matthews", title = "Illinois")
Plots.plot!(-demean_standardize(permed_cis[2, :]), label = "Pseudo-Structural")
cor(coincident[:, 2], permed_cis[2, :])

# Indiana
Plots.plot(demean_standardize(coincident[:, 3]), label = "Crone and Clayton-Matthews")
Plots.plot!(-demean_standardize(permed_cis[3, :]), label = "Pseudo-Structural", title = "Indiana")
cor(coincident[:, 3], permed_cis[3, :])

# Minnesota
Plots.plot(demean_standardize(coincident[:, 5]), label = "Crone and Clayton-Matthews")
Plots.plot!(-demean_standardize(permed_cis[4, :]), label = "Pseudo-Structural", title = "Minnesota")
cor(coincident[:, 5], permed_cis[4, :])

# North Dakota
Plots.plot(demean_standardize(coincident[:, 6]), label = "Crone and Clayton-Matthews")
Plots.plot!(-demean_standardize(permed_cis[5, :]), label = "Pseudo-Structural", title = "North Dakota")
cor(coincident[:, 6], permed_cis[5, :])
# Plots.plot!(demean_standardize(employment[:, 5]), label = "Employment")
# ND experienced an oil boom between 2006, peaked in 2012, and crashed in 2014.
# can plot it as
# Plots.plot(ut_employment[:, 5])

# Ohio
Plots.plot(demean_standardize(coincident[:, 7]), label = "Crone and Clayton-Matthews")
Plots.plot!(-demean_standardize(permed_cis[6, :]), label = "Pseudo-Structural", title = "Ohio")
cor(coincident[:, 7], permed_cis[6, :])

# South Dakota
Plots.plot(demean_standardize(coincident[:, 8]), label = "Crone and Clayton-Matthews")
Plots.plot!(-demean_standardize(permed_cis[7, :]), label = "Pseudo-Structural", title = "South Dakota")
cor(coincident[:, 8], permed_cis[7, :])

# Delta method using the asymptotic variance for delta and gamma
dg_var = inv(res.hess_est)[1:13, 1:13]

# I need delta star and gamma star
delta_star = res.delta_est[6:end, :]
gamma_star = res.gamma_est[end, :]'

# I want the asymptotic variance of kron(gamma_star, delta_star)'
m, n = size(reshape(gamma_star', :, 1))
p, q = size(delta_star')
commutation_mat = comm_mat(q, m)
large_comm = kron(kron(I(n), commutation_mat), I(p))
partial_gamma = kron(I(3), vec(delta_star'))
partial_delta = kron(vec(gamma_star'), I(10))
J = large_comm * hcat(partial_delta, partial_gamma)

# Delta result
delta_method_variance = J * dg_var * J'
delta_method_stderr = sqrt.(diag(delta_method_variance))

theta = vec(kron(gamma_star', delta_star'))
tstats = theta ./ delta_method_stderr

pvals = 2 .* (1 .- cdf.(Normal(), abs.(tstats)))

sig_5 = abs.(tstats) .> 1.372
hcat(theta, delta_method_stderr, tstats, sig_5)
# It is really only the second variable in the list that is significant

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









