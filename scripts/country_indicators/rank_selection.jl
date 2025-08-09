using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250723)
using RCall

R"""
source("r_helpers.R")
"""

matdata = load(datadir("./country_indicators/globaldata.jld2"), "matdata");
vecdata = tenmat(matdata, row=[1, 2])
pdata = permutedims(matdata, (3, 1, 2))

dimvals = collect(size(matdata[:, :, 1]))
cen_data = vecdata .- mean(vecdata, dims=2)

smallicest = rank_selection(cen_data, dimvals; iters=1000, pmax=4)
# AIC selects 3,5,1
# BIC selects 2,5,1
# HQC selects 3,5,1

# AIC selects 2,3,1
# BIC selects 2,1,1
# HQC selects 2,2,1


small_bench = R"""
d1 = $dimvals[1]
d2 = $dimvals[2]
data <- $pdata
selected_rank <- r_rank_selection(data, d1, d2)
"""
@rget selected_rank

tmp = R"""
est <- matAR.RR.est(data, method = "RRMLE", k1 = 3, k2 = 1)
sd <- matAR.RR.se(est$A1, est$A2, 3, 1, method = "RRMLE", Sigma1=est$Sig1, Sigma2=est$Sig2)
"""
@rget est
@rget sd



omega = res18.omega
sigma1 = res18.sigma1_est
sigma2 = res18.sigma2_est
det(sigma1)
det(sigma2)
kron_term = kron(sigma2, sigma1)
det(kron_term)
eigvals(kron_term)
mid_term = inv(omega * kron(sigma2, sigma1) * omega')
eigvals(mid_term)
perm_mat = perm_matrix(dimvals, [2,5])

std(cen_data, dims = 2)

function system_parameters(dimvals, ranks; p=1)
    first_term = ranks[1] * (dimvals[1] - ranks[1])
    second_term = ranks[2] * (dimvals[2] - ranks[2])
    third_term = p * ranks[1] * dimvals[1] + p * ranks[2] * dimvals[2] - p
    num_ll1 = Int(dimvals[1] * (dimvals[1] + 1) / 2) - 1
    num_ll2 = Int(dimvals[2] * (dimvals[2] + 1) / 2)

    return first_term + second_term + num_ll1 + num_ll2 + third_term - 1
end



aic(ll::Real, numpars::Int) = -2 * ll + (2 * numpars)
bic(ll::Real, numpars::Int, obs::Int) = -2 * ll + (numpars2 * log(obs))
hqc(ll::Real, numpars::Int, obs::Int) = -2 * ll + (numpars * 2 * log(log(obs)))
ll1 = -res3.res.minimum
ll2 = -res20.res.minimum
ll3 = -res1.res.minimum
numpars1 = system_parameters(dimvals, [3,1])
numpars2 = system_parameters(dimvals, [4,5])
numpars3 = system_parameters(dimvals, [4,5]; p=2)
obs = 99
bic1 = bic(ll1, numpars1, obs)
bic2 = bic(ll2, numpars2, obs)
bic3 = bic(ll3, numpars2, obs)

numpars3 = system_parameters([3,4], [3,4]; p=3)



res1 = comovement_reg(cen_data, dimvals, [1, 1]; iters=1000, p=4)
res2 = comovement_reg(cen_data, dimvals, [2, 1]; iters=1000)
res3 = comovement_reg(cen_data, dimvals, [3, 1]; iters=1000)
res4 = comovement_reg(cen_data, dimvals, [4, 1]; iters=1000)
res5 = comovement_reg(cen_data, dimvals, [1, 2]; iters=1000)
res6 = comovement_reg(cen_data, dimvals, [2, 2]; iters=1000)
res7 = comovement_reg(cen_data, dimvals, [3, 2]; iters=1000)
res8 = comovement_reg(cen_data, dimvals, [4, 2]; iters=1000)
res9 = comovement_reg(cen_data, dimvals, [1, 3]; iters=1000)
res10 = comovement_reg(cen_data, dimvals, [2, 3]; iters=1000)
res11 = comovement_reg(cen_data, dimvals, [3, 3]; iters=1000)
res12 = comovement_reg(cen_data, dimvals, [4, 3]; iters=1000)
res13 = comovement_reg(cen_data, dimvals, [1, 4]; iters=1000)
res14 = comovement_reg(cen_data, dimvals, [2, 4]; iters=1000)
res15 = comovement_reg(cen_data, dimvals, [3, 4]; iters=1000, p=4)
res16 = comovement_reg(cen_data, dimvals, [4, 4]; iters=1000)
res17 = comovement_reg(cen_data, dimvals, [1, 5]; iters=1000)
res18 = comovement_reg(cen_data, dimvals, [2, 5]; iters=1000)
res19 = comovement_reg(cen_data, dimvals, [3, 5]; iters=1000)
res20 = comovement_reg(cen_data, dimvals, [4, 5]; iters=1000, p=4)
res.res.minimum

res.delta_est
res.delta_stderr
res.gamma_est
res.gamma_stderr

s1 = res2.sigma1_est ./ norm(res2.sigma1_est)
s2 = res2.sigma2_est .* norm(res2.sigma1_est)

res2.delta_est
res2.delta_stderr

res2.gamma_est
res2.gamma_stderr
