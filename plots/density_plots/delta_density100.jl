using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/delta_cov_results100.jld2"))

correct_delta = density_results["correct_delta"]
under_delta = density_results["under_delta"]
over_delta = density_results["over_delta"]
delta_true = density_results["delta_true"]

h1 = density(correct_delta[1, :], legend=false, ylabel="Density")
density!(under_delta[1, :])
density!(over_delta[1, :])
vline!([delta_true[2]])
savefig(h1, "delta_density100_1.png")

h2 = density(correct_delta[2, :], legend=false)
density!(under_delta[2, :])
density!(over_delta[2, :])
vline!([delta_true[3]])
savefig(h2, "delta_density100_2.png")

h3 = density(correct_delta[3, :], label="Correct rank", legendfontsize=10)
density!(under_delta[3, :], label="Underestimated rank")
density!(over_delta[3, :], label="Overestimated rank")
vline!([delta_true[4]], label="True value")
savefig(h3, "delta_density100_3.png")
