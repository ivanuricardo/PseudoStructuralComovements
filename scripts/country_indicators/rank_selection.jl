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

smallicest = rank_selection(vecdata, dimvals; iters=1000)

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

res1 = comovement_reg(cen_data, dimvals, [1, 1]; iters=1000)
res2 = comovement_reg(cen_data, dimvals, [2, 1]; iters=1000)
res3new = comovement_reg(cen_data, dimvals, [3, 1]; iters=1000)
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
res15 = comovement_reg(cen_data, dimvals, [3, 4]; iters=1000)
res16 = comovement_reg(cen_data, dimvals, [4, 4]; iters=1000)
res17 = comovement_reg(cen_data, dimvals, [1, 5]; iters=1000)
res18 = comovement_reg(cen_data, dimvals, [2, 5]; iters=1000)
res19 = comovement_reg(cen_data, dimvals, [3, 5]; iters=1000)
res20 = comovement_reg(cen_data, dimvals, [4, 5]; iters=1000)
res.res.minimum

res.delta_est
res.delta_stderr
res.gamma_est
res.gamma_stderr


