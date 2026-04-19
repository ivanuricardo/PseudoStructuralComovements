using DrWatson
@quickactivate :PseudoStructuralComovements

density_results = load(datadir("coverage/34/gamma_comparison_results100.jld2"))
correct_gamma = density_results["correct_gamma"]
under_gamma   = density_results["under_gamma"]
over_gamma    = density_results["over_gamma"]
correct_rrmar = density_results["correct_rrmar"]
under_rrmar   = density_results["under_rrmar"]
over_rrmar    = density_results["over_rrmar"]
gamma_true    = density_results["gamma_true"]

correct_cov = density_results["correct_cov"]
under_cov = density_results["under_cov"]
over_cov = density_results["over_cov"]

h1 = StatsPlots.density(
    correct_gamma[1, :];
    legend=false,
    ylabel="Density",
    linewidth=3,
    yguidefont=16   # <-- increase y-axis label font size
)
StatsPlots.density!(under_gamma[1, :]; linewidth=3)
StatsPlots.density!(over_gamma[1, :]; linewidth=3)
vline!([gamma_true[2]]; linewidth=3)
savefig(h1, "plots/density_plots/34/gamma_density100_1.pdf")

h2 = StatsPlots.density(
    correct_gamma[2, :];
    legend=false,
    linewidth=3,
)
StatsPlots.density!(under_gamma[2, :]; linewidth=3)
StatsPlots.density!(over_gamma[2, :]; linewidth=3)
vline!([gamma_true[3]]; linewidth=3)
savefig(h2, "plots/density_plots/34/gamma_density100_2.pdf")

h3 = StatsPlots.density(
    correct_gamma[3, :];
    linewidth=3,
    label="Correct",
    legendfontsize=13,
)
StatsPlots.density!(under_gamma[3, :]; label="Underestimated", linewidth=3)
StatsPlots.density!(over_gamma[3, :]; label="Overestimated", linewidth=3)
vline!([gamma_true[4]]; linewidth=3,label="True value")
savefig(h3, "plots/density_plots/34/gamma_density100_3.pdf")

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

r1 = StatsPlots.density(
    correct_rrmar[1, :];
    legend=false,
    ylabel="Density",
    linewidth=3,
    yguidefont=16   # <-- increase y-axis label font size
)
StatsPlots.density!(under_rrmar[1, :]; linewidth=3)
StatsPlots.density!(over_rrmar[1, :]; linewidth=3)
vline!([gamma_true[2]]; linewidth=3)
savefig(r1, "plots/density_plots/34/rrmar_gamma_density100_1.pdf")

r2 = StatsPlots.density(
    correct_rrmar[2, :];
    legend=false,
    linewidth=3,
    label="Correct",
    legendfontsize=13,
)
StatsPlots.density!(under_rrmar[2, :]; label="Underestimated", linewidth=3)
StatsPlots.density!(over_rrmar[2, :]; label="Overestimated", linewidth=3)
vline!([gamma_true[3]]; linewidth=3,label="True value")
savefig(r2, "plots/density_plots/34/rrmar_gamma_density100_2.pdf")

r3 = StatsPlots.density(
    correct_rrmar[3, :];
    linewidth=3,
    label="Correct",
    legendfontsize=13,
)
StatsPlots.density!(under_rrmar[3, :]; label="Underestimated", linewidth=3)
StatsPlots.density!(over_rrmar[3, :]; label="Overestimated", linewidth=3)
vline!([gamma_true[4]]; linewidth=3,label="True value")
savefig(r3, "plots/density_plots/34/rrmar_gamma_density100_3.pdf")
