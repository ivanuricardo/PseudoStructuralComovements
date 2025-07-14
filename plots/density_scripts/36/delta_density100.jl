using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/36/delta_cov_results100.jld2"))

correct_delta = density_results["correct_delta"]
under_delta = density_results["under_delta"]
over_delta = density_results["over_delta"]
delta_true = density_results["delta_true"]

h1 = density(correct_delta[1, :], legend=false, ylabel="Density")
density!(under_delta[1, :])
density!(over_delta[1, :])
vline!([delta_true[2]])
savefig(h1, "../density_plots/36/delta_density100_1.png")

h2 = density(correct_delta[2, :], legend=false)
density!(under_delta[2, :])
density!(over_delta[2, :])
vline!([delta_true[3]])
savefig(h2, "../density_plots/36/delta_density100_2.png")

