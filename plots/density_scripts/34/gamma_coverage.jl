using DrWatson
@quickactivate :PseudoStructuralComovements

density_results100 = load(datadir("coverage/34/gamma_cov_results100.jld2"))
density_results250 = load(datadir("coverage/34/gamma_cov_results250.jld2"))

correct_cov100 = density_results100["correct_cov"]
under_cov100 = density_results100["under_cov"]
over_cov100 = density_results100["over_cov"]
correct100 = vec(sum(correct_cov100, dims=2))
under100 = vec(sum(under_cov100, dims=2))
over100 = vec(sum(over_cov100, dims=2))

correct_cov250 = density_results250["correct_cov"]
under_cov250 = density_results250["under_cov"]
over_cov250 = density_results250["over_cov"]
correct250 = vec(sum(correct_cov250, dims=2))
under250 = vec(sum(under_cov250, dims=2))
over250 = vec(sum(over_cov250, dims=2))

# Labels
gammas = ["γ₁", "γ₂", "γ₃"]
delta_ranks = ["Correct", "Overestimated", "Underestimated"]
sample_sizes = ["T = 100", "T = 250"]

# Coverage results (in percentages)
# Format: gamma_estimate[delta_rank][sample_size]

# Now build coverage_data by looping over the three γ’s:
coverage_data = Dict{String,Vector{Vector{Float64}}}()
for (i, γ) in enumerate(gammas)
    coverage_data[γ] = [
        [correct100[i], over100[i], under100[i]],  # T = 100
        [correct250[i], over250[i], under250[i]]   # T = 250
    ]
end

# Convert to percentage
coverage_data = Dict(k => [map(x -> x / 10, row) for row in v] for (k, v) in coverage_data)

# Prepare plotting data
bar_width = 0.25
x_positions = [1, 2, 3]  # Delta rank groups
offsets = [-bar_width, 0.0, bar_width]  # For each gamma estimator

fig = Figure(
    size=(800, 500),
    backgroundcolor=:transparent     # make the figure background invisible
)

for (j, T) in enumerate(sample_sizes)
    ax = Axis(fig[j, 1],
        title=T,
        ylabel="Coverage (%)",
        xticks=(x_positions, delta_ranks),
        limits=(nothing, (90, 100)),
        yticks=90:2:100,
        backgroundcolor=:transparent,
    )

    for (i, est) in enumerate(gammas)
        coverage = coverage_data[est][j]
        barpositions = x_positions .+ offsets[i]
        barplot!(barpositions, coverage;
            width=bar_width,
            label=est,
        )
    end

    axislegend(ax; position=:rt, backgroundcolor=:transparent)
end

save("../../density_plots/34/gamma_coverage_plot.png", fig)
