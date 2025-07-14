using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/34/gamma_cov_results250.jld2"))

correct_gamma = density_results["correct_gamma"]
under_gamma = density_results["under_gamma"]
over_gamma = density_results["over_gamma"]
gamma_true = density_results["gamma_true"]

h1 = density(correct_gamma[1, :], legend=false, ylabel="Density")
density!(under_gamma[1, :])
density!(over_gamma[1, :])
vline!([gamma_true[2]])
savefig(h1, "../density_plots/34/gamma_density250_1.png")

h2 = density(correct_gamma[2, :], legend=false)
density!(under_gamma[2, :])
density!(over_gamma[2, :])
vline!([gamma_true[3]])
savefig(h2, "../density_plots/34/gamma_density250_2.png")

h3 = density(correct_gamma[3, :], legend=false)
density!(under_gamma[3, :], label="Underestimated rank")
density!(over_gamma[3, :], label="Overestimated rank")
vline!([gamma_true[4]], label="True value")
savefig(h3, "../density_plots/34/gamma_density250_3.png")
