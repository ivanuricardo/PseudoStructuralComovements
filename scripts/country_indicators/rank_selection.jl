using DrWatson
@quickactivate :PseudoStructuralComovements

matdata = load(datadir("./country_indicators/globaldata.jld2"), "matdata");
vecdata = tenmat(matdata, row=[1, 2])

dimvals = collect(size(matdata[:, :, 1]))

smallicest = rank_selection(vecdata, dimvals; iters=1000)
