using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250723)

matdata = load(datadir("./country_indicators/globaldata.jld2"), "matdata");
vecdata = tenmat(matdata, row=[1, 2])

dimvals = collect(size(matdata[:, :, 1]))
cen_data = vecdata .- mean(vecdata, dims=2)

smallicest = rank_selection(vecdata, dimvals; iters=1000)

res = comovement_reg(cen_data, dimvals, [3, 1]; iters=1000)
res.res.minimum

res.delta_est
res.delta_stderr
res.gamma_est
res.gamma_stderr


