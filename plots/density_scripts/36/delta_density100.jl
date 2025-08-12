using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/36/delta_cov_results100.jld2"))

correct_delta = density_results["correct_delta"]
under_delta = density_results["under_delta"]
delta_true = density_results["delta_true"]

h1 = Plots.density(correct_delta[1, :], legend=false, ylabel="Density")
Plots.density!(under_delta[1, :])
Plots.vline!([delta_true[2]])
savefig(h1, "../../density_plots/36/delta_density100_1_36.png")

h2 = Plots.density(correct_delta[2, :], label="Correct gamma", legendfontsize=9)
Plots.density!(under_delta[2, :], label="Underestimated delta")
Plots.vline!([delta_true[3]], label = "True Value")
savefig(h2, "../../density_plots/36/delta_density100_2_36.png")

