using DrWatson
@quickactivate :PseudoStructuralComovements
using CSV, Tables, SparseArrays, Distributions, XLSX
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

state_names = ["IA", "IL", "IN", "MN", "ND", "OH", "SD"]
rawdata = XLSX.readdata(datadir("./state_indexes/reguib_northcentral.xlsx"), "Sheet1!A2:S459")
vecdata = Float64.(rawdata[:, 2:end])'
coincident = vecdata[1:2:end, 230:end-2]'

start_date = 133  # Jan 2001
covid_start = 12  # Dec 2019

ut_employment = load_series("employment")[start_date:(end-covid_start), :]
ut_unemployment = load_series("unemployment")[start_date:(end-covid_start), :]
ut_hours = load_series("hours")[1:(end-covid_start), :]  # already starts in Jan

employment = transform(ut_employment)
unemployment = transform(ut_unemployment; type = "diff")
hours = transform(ut_hours)

file = XLSX.readxlsx(datadir("updated_states/wages.xlsx"))
sheet = file["Table"]
ut_wages = Float64.(XLSX.getdata(sheet, "O7:CL13"))'
# Rearrange to be alphabetical and align with other series
# iowa, illinois, indiana, minnesota, north dakota, ohio, south dakota
rearranged_wages = ut_wages[:, [3,1,2,4,5,6,7]]
monthly_wages = denton(rearranged_wages, ut_employment)
wages = transform(monthly_wages)

catted_data = cat(unemployment', hours', wages', employment'; dims = 3)
tendata = permutedims(catted_data, (1,3,2))
perm_states = [7,4,5,1,6,2,3]
# SD, MN, ND, IA, OH, IL, IN
# perm_states = [1,2,5,6,7,4,3]
# IA, IL, MN, OH, SD, ND, IN
rearranged_tendata = tendata[perm_states, :, :]
matdata = vectorize(rearranged_tendata)

n1, n2, obs = size(tendata)
dimvals = [n1, n2]
# res = comovement_reg(matdata, dimvals, [2, 1]; iters=1000, p=1)
# save the results
# save(datadir("updated_states/coincident_results.jld2"), Dict("res" => res))
# load the results
loaded_results = load(datadir("updated_states/coincident_results.jld2"))
res = loaded_results["res"]

# julia> res.delta_est
# 7×5 Matrix{Float64}:
#   1.0              0.0              0.0              0.0              0.0
#   0.0              1.0              0.0              0.0              0.0
#   0.0              0.0              1.0              0.0              0.0
#   0.0              0.0              0.0              1.0              0.0
#   0.0              0.0              0.0              0.0              1.0
#   0.110(0.101)     0.022(0.096)    -0.084(0.113)     0.074(0.080)     0.036(0.085)
#  -0.475(0.094)    -0.611(0.090)    -0.198(0.105)    -0.635(0.075)    -0.931(0.079)

# South Dakota co-moves Indiana, not so much with Illinois
# Minnesota co-moves with Indiana, not so much with Illinois
# North Dakota doesn't co-move with any
# Iowa co-moves with Indiana
# Ohio co-moves with Indiana
# So there seems to be a grouping of countries with
# Ohio, Indiana, Iowa, Minnesota, South Dakota
# North Dakota and Illinois don't seem to co-move with anything


# julia> res.gamma_est
# 4×3 Matrix{Float64}:
#   1.0              0.0               0.0
#   0.0              1.0               0.0
#   0.0              0.0               1.0
#  -0.599(0.057)     58.619(6.562)    -0.258(0.412)

# employment co-moves with wages
# unemployment co-moves with wages
# hours co-moves with wages
# All co-move significantly, which makes sense and can create the coincident
# indicator


################################################################################
# Plots
################################################################################

# cis order is
# SD, MN, ND, IA, OH, IL, IN
# coincident order is 
# IA, IL, IN, MI, MN, ND, OH, SD, WI
# we omit MI and WI
factors = kron(res.u4_est, res.u3_est)' * matdata
u1 = nullspace(res.delta_est') * inv(nullspace(res.delta_est')[6:end, :])
cis = u1 * factors
permed_cis = cis[invperm(perm_states), :]

# Iowa
Plots.plot(demean_standardize(coincident[:, 1]), label = "Crone and Clayton-Matthews", title = "Iowa")
Plots.plot!(demean_standardize(permed_cis[1, :]), label = "Pseudo-Structural")

# Illinois
Plots.plot(demean_standardize(coincident[:, 2]), label = "Crone and Clayton-Matthews", title = "Illinois")
Plots.plot!(demean_standardize(permed_cis[2, :]), label = "Pseudo-Structural")

# Indiana
Plots.plot(demean_standardize(coincident[:, 3]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[3, :]), label = "Pseudo-Structural", title = "Indiana")

# Minnesota
Plots.plot(demean_standardize(coincident[:, 5]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[4, :]), label = "Pseudo-Structural", title = "Minnesota")

# North Dakota
Plots.plot(demean_standardize(coincident[:, 6]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[5, :]), label = "Pseudo-Structural", title = "North Dakota")
# Plots.plot!(demean_standardize(employment[:, 5]), label = "Employment")
# ND experienced an oil boom between 2006, peaked in 2012, and crashed in 2014.
# can plot it as
# Plots.plot(ut_employment[:, 5])

# Ohio
Plots.plot(demean_standardize(coincident[:, 7]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[6, :]), label = "Pseudo-Structural", title = "Ohio")

# South Dakota
Plots.plot(demean_standardize(coincident[:, 8]), label = "Crone and Clayton-Matthews")
Plots.plot!(demean_standardize(permed_cis[7, :]), label = "Pseudo-Structural", title = "South Dakota")

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












