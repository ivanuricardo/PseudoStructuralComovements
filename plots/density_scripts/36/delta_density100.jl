using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/36/delta_cov_results100.jld2"))

correct_delta = density_results["correct_delta"]
under_delta = density_results["under_delta"]
delta_true = density_results["delta_true"]

h1 = StatsPlots.density(
    correct_delta[1, :];
    legend=false,
    ylabel="Density",
    linewidth=3,
    yguidefont=16   # <-- increase y-axis label font size
)
StatsPlots.density!(under_delta[1, :]; linewidth=3)
Plots.vline!([delta_true[2]], linewidth=3)
savefig(h1, "../../density_plots/36/delta_density100_1_36.png")

h2 = StatsPlots.density(
    correct_delta[2, :];
    linewidth=3,
    label="Correct",
    legendfontsize=13,
)
StatsPlots.density!(under_delta[2, :]; label="Underestimated", linewidth=3)
vline!([delta_true[3]]; linewidth=3,label="True value")
savefig(h2, "../../density_plots/36/delta_density100_2_36.png")

