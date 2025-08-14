using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250723)
using RCall, XLSX

R"""
source("r_helpers.R")
"""

matdata = XLSX.readdata(datadir("./state_indexes/reguib_northcentral.xlsx"), "Sheet1!A2:S459")
vecdata = Float64.(matdata[:, 2:end])'

dimvals = [2,9]
cen_data = vecdata .- mean(vecdata, dims=2)

icest = rank_selection(cen_data, dimvals; iters=1000, pmax=3)

# AIC selects ranks 2,9 with 3 lags
# BIC and HQC selects ranks 2,9 with 2 lags


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

hess_est = res1.hess_est



aic(ll::Real, numpars::Int) = -2 * ll + (2 * numpars)
bic(ll::Real, numpars::Int, obs::Int) = -2 * ll + (numpars2 * log(obs))
hqc(ll::Real, numpars::Int, obs::Int) = -2 * ll + (numpars * 2 * log(log(obs)))

ll1 = -res1.res.minimum
ll2 = -res2.res.minimum
ll3 = -res3.res.minimum
ll4 = -res4.res.minimum
ll5 = -res5.res.minimum
ll6 = -res6.res.minimum
ll7 = -res7.res.minimum
ll8 = -res8.res.minimum
ll9 = -res9.res.minimum
ll10 = -res10.res.minimum
ll11 = -res11.res.minimum
ll12 = -res12.res.minimum
ll13 = -res13.res.minimum
ll14 = -res14.res.minimum
ll15 = -res15.res.minimum
ll16 = -res16.res.minimum
ll17 = -res17.res.minimum
ll18 = -res18.res.minimum


numpars1 = system_parameters(dimvals, [1, 1]; p=3)
numpars2 = system_parameters(dimvals, [2, 1]; p=3)
numpars3 = system_parameters(dimvals, [1, 2]; p=3)
numpars4 = system_parameters(dimvals, [2, 2]; p=3)
numpars5 = system_parameters(dimvals, [1, 3]; p=3)
numpars6 = system_parameters(dimvals, [2, 3]; p=2)
numpars7 = system_parameters(dimvals, [1, 4]; p=2)
numpars8 = system_parameters(dimvals, [2, 4]; p=2)
numpars9 = system_parameters(dimvals, [1, 5]; p=2)
numpars10 = system_parameters(dimvals, [2,5]; p=2)
numpars11 = system_parameters(dimvals, [1,6]; p=2)
numpars12 = system_parameters(dimvals, [2,6]; p=2)
numpars13 = system_parameters(dimvals, [1,7]; p=2)
numpars14 = system_parameters(dimvals, [2,7]; p=2)
numpars15 = system_parameters(dimvals, [1,8]; p=2)
numpars16 = system_parameters(dimvals, [2,8]; p=2)
numpars17 = system_parameters(dimvals, [1,9]; p=2)
numpars18 = system_parameters(dimvals, [2,9]; p=2)

obs = 457
bic1 = bic(ll1, numpars1, obs)
bic2 = bic(ll2, numpars2, obs)
bic3 = bic(ll3, numpars3, obs)
bic4 = bic(ll4, numpars4, obs)
bic5 = bic(ll5, numpars5, obs)
bic6 = bic(ll6, numpars6, obs)
bic7 = bic(ll7, numpars7, obs)
bic8 = bic(ll8, numpars8, obs)
bic9 = bic(ll9, numpars9, obs)
bic10 = bic(ll10, numpars10, obs)
bic11 = bic(ll11, numpars11, obs)
bic12 = bic(ll12, numpars12, obs)
bic13 = bic(ll13, numpars13, obs)
bic14 = bic(ll14, numpars14, obs)
bic15 = bic(ll15, numpars15, obs)
bic16 = bic(ll16, numpars16, obs)
bic17 = bic(ll17, numpars17, obs)
bic18 = bic(ll18, numpars18, obs)

# BIC chooses ranks 2,9 with bic value -14252.4629142954 (ll value is 7527.39)

res1 = comovement_reg(cen_data, dimvals, [1, 1]; iters=1000, p=3)
res2 = comovement_reg(cen_data, dimvals, [2, 1]; iters=1000, p=3)
res3 = comovement_reg(cen_data, dimvals, [1, 2]; iters=1000, p=3)
res4 = comovement_reg(cen_data, dimvals, [2, 2]; iters=1000, p=2)
res5 = comovement_reg(cen_data, dimvals, [1, 3]; iters=1000, p=2)
res6 = comovement_reg(cen_data, dimvals, [2, 3]; iters=1000, p=2)
res7 = comovement_reg(cen_data, dimvals, [1, 4]; iters=1000, p=2)
res8 = comovement_reg(cen_data, dimvals, [2, 4]; iters=1000, p=2)
res9 = comovement_reg(cen_data, dimvals, [1, 5]; iters=1000, p=2)

# From res10, res8 is the best so far with -7.217803e+03
res10 = comovement_reg(cen_data, dimvals, [2, 5]; iters=1000, p=2)
res11 = comovement_reg(cen_data, dimvals, [1, 6]; iters=1000, p=2)
res12 = comovement_reg(cen_data, dimvals, [2, 6]; iters=1000, p=2)
res13 = comovement_reg(cen_data, dimvals, [1, 7]; iters=1000, p=2)
res14 = comovement_reg(cen_data, dimvals, [2, 7]; iters=1000, p=2)
res15 = comovement_reg(cen_data, dimvals, [1, 8]; iters=1000, p=2)
res16 = comovement_reg(cen_data, dimvals, [2, 8]; iters=1000, p=2)
res17 = comovement_reg(cen_data, dimvals, [1, 9]; iters=1000, p=2)
res18 = comovement_reg(cen_data, dimvals, [2, 9]; iters=1000, p=2)

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
