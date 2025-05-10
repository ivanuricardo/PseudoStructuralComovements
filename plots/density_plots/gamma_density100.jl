using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/gamma_cov_results100.jld2"))

correct_gamma = density_results["correct_gamma"]
under_gamma = density_results["under_gamma"]
over_gamma = density_results["over_gamma"]
gamma_true = density_results["gamma_true"]

sum(density_results["correct_cov"], dims=2)
sum(density_results["under_cov"], dims=2)
sum(density_results["over_cov"], dims=2)

h1 = density(correct_gamma[1, :], label="Correct rank")
density!(under_gamma[1, :], label="Underestimated rank")
density!(over_gamma[1, :], label="Overestimated rank")
vline!([gamma_true[2]], label="True value")
vline!([mean(correct_gamma[1, :])], label="Correct rank mean")
vline!([mean(under_gamma[1, :])], label="Under rank mean")
vline!([mean(over_gamma[1, :])], label="Over rank mean")

h2 = density(correct_gamma[2, :], label="Correct rank")
density!(under_gamma[2, :], label="Underestimated rank")
density!(over_gamma[2, :], label="Overestimated rank")
vline!([gamma_true[3]], label="True value")

h3 = density(correct_gamma[3, :], label="Correct rank")
density!(under_gamma[3, :], label="Underestimated rank")
density!(over_gamma[3, :], label="Overestimated rank")
vline!([gamma_true[4]], label="True value")
