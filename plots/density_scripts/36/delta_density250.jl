using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/36/delta_cov_results250.jld2"))

correct_delta = density_results["correct_delta"]
under_delta = density_results["under_delta"]
delta_true = density_results["delta_true"]

h1 = Plots.density(correct_delta[1, :], legendfontsize=9, ylabel="Density", label = "Correct gamma")
Plots.density!(under_delta[1, :], label = "Underestimated gamma")
Plots.vline!([delta_true[2]], label = "True value")
savefig(h1, "../../density_plots/36/delta_density250_1.png")

h2 = Plots.density(correct_delta[2, :], legend=false, label = "Correct gamma")
Plots.density!(under_delta[2, :], label = "Underestimated gamma")
Plots.vline!([delta_true[3]], label = "True Value")
savefig(h2, "../../density_plots/36/delta_density250_2.png")
