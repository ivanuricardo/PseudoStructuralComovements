using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/34/delta_cov_results100.jld2"))

correct_delta = density_results["correct_delta"]
under_delta = density_results["under_delta"]
over_delta = density_results["over_delta"]
delta_true = density_results["delta_true"]

correct_cov = density_results["correct_cov"]
under_cov = density_results["under_cov"]
over_cov = density_results["over_cov"]

h1 = StatsPlots.density(
    correct_delta[1, :];
    legend=false,
    ylabel="Density",
    linewidth=3,
    yguidefont=16   # <-- increase y-axis label font size
)
StatsPlots.density!(under_delta[1, :]; linewidth=3)
StatsPlots.density!(over_delta[1, :]; linewidth=3)
vline!([delta_true[2]]; linewidth=3)
savefig(h1, "../../density_plots/34/delta_density100_1.png")

h2 = StatsPlots.density(
    correct_delta[2, :];
    linewidth=3,
    label="Correct",
    legendfontsize=13,
)
StatsPlots.density!(under_delta[2, :]; label="Underestimated", linewidth=3)
StatsPlots.density!(over_delta[2, :]; label="Overestimated", linewidth=3)
vline!([delta_true[3]]; linewidth=3,label="True value")
savefig(h2, "../../density_plots/34/delta_density100_2.png")

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
