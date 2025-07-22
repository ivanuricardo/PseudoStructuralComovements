using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/34/gamma_cov_results100.jld2"))

correct_gamma = density_results["correct_gamma"]
under_gamma = density_results["under_gamma"]
over_gamma = density_results["over_gamma"]
gamma_true = density_results["gamma_true"]

correct_cov = density_results["correct_cov"]
under_cov = density_results["under_cov"]
over_cov = density_results["over_cov"]

h1 = density(correct_gamma[1, :], legend=false, ylabel="Density")
density!(under_gamma[1, :])
density!(over_gamma[1, :])
vline!([gamma_true[2]])
savefig(h1, "../../density_plots/34/gamma_density100_1.png")

h2 = density(correct_gamma[2, :], legend=false)
density!(under_gamma[2, :])
density!(over_gamma[2, :])
vline!([gamma_true[3]])
savefig(h2, "../../density_plots/34/gamma_density100_2.png")

h3 = density(correct_gamma[3, :], label="Correct delta", legendfontsize=9)
density!(under_gamma[3, :], label="Underestimated delta")
density!(over_gamma[3, :], label="Overestimated delta")
vline!([gamma_true[4]], label="True value")
savefig(h3, "../../density_plots/34/gamma_density100_3.png")

correct_coverages = sum(correct_cov, dims=2)
over_coverages = sum(over_cov, dims=2)
under_coverages = sum(under_cov, dims=2)

println("First gamma estimate coverages:")
println("given correct delta rank: $(correct_coverages[1])")
println("given overestimated delta rank: $(over_coverages[1])")
println("given underestimated delta rank: $(under_coverages[1])")
println("_________________________________________")

println("Second gamma estimate coverages:")
println("given correct delta rank: $(correct_coverages[2])")
println("given overestimated delta rank: $(over_coverages[2])")
println("given underestimated delta rank: $(under_coverages[2])")
println("_________________________________________")

println("Third gamma estimate coverages:")
println("given correct delta rank: $(correct_coverages[3])")
println("given overestimated delta rank: $(over_coverages[3])")
println("given underestimated delta rank: $(under_coverages[3])")
