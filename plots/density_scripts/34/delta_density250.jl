using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/34/delta_cov_results250.jld2"))

correct_delta = density_results["correct_delta"]
under_delta = density_results["under_delta"]
over_delta = density_results["over_delta"]
delta_true = density_results["delta_true"]

correct_cov = density_results["correct_cov"]
under_cov = density_results["under_cov"]
over_cov = density_results["over_cov"]

h1 = StatsPlots.density(correct_delta[1, :], legend=false, ylabel="Density")
StatsPlots.density!(under_delta[1, :])
StatsPlots.density!(over_delta[1, :])
vline!([delta_true[2]])
savefig(h1, "../../density_plots/34/delta_density250_1.png")

h2 = StatsPlots.density(correct_delta[2, :], legend=false)
StatsPlots.density!(under_delta[2, :])
StatsPlots.density!(over_delta[2, :])
vline!([delta_true[3]])
savefig(h2, "../../density_plots/34/delta_density250_2.png")

correct_coverages = sum(correct_cov, dims=2)
over_coverages = sum(over_cov, dims=2)
under_coverages = sum(under_cov, dims=2)

println("First delta estimate coverages:")
println("given correct gamma rank: $(correct_coverages[1])")
println("given overestimated gamma rank: $(over_coverages[1])")
println("given underestimated gamma rank: $(under_coverages[1])")
println("_________________________________________")

println("Second delta estimate coverages:")
println("given correct gamma rank: $(correct_coverages[2])")
println("given overestimated gamma rank: $(over_coverages[2])")
println("given underestimated gamma rank: $(under_coverages[2])")
