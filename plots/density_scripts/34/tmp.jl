using DrWatson
@quickactivate :PseudoStructuralComovements

# Labels
estimators = ["γ₁", "γ₂", "γ₃"]
delta_ranks = ["Correct", "Overestimated", "Underestimated"]
sample_sizes = ["T = 100", "T = 250"]

# Coverage results (in percentages)
# Format: gamma_estimate[delta_rank][sample_size]
coverage_data = Dict(
    "γ₁" => [  # T=100, T=250
        [942, 952, 958],  # T=100
        [949, 952, 934]   # T=250
    ],
    "γ₂" => [
        [942, 944, 954],
        [943, 943, 950]
    ],
    "γ₃" => [
        [938, 938, 944],
        [952, 950, 957]
    ]
)

# Convert to percentage
coverage_data = Dict(k => [map(x -> x / 10, row) for row in v] for (k, v) in coverage_data)

# Prepare plotting data
bar_width = 0.25
x_positions = [1, 2, 3]  # Delta rank groups
offsets = [-bar_width, 0.0, bar_width]  # For each gamma estimator

fig = Figure(resolution=(800, 500))

for (j, T) in enumerate(sample_sizes)
    ax = Axis(fig[j, 1],
        title=T,
        ylabel="Coverage (%)",
        xticks=(x_positions, delta_ranks),
        limits=(nothing, (90, 100)),
        yticks=90:2:100
    )

    for (i, est) in enumerate(estimators)
        coverage = coverage_data[est][j]
        barpositions = x_positions .+ offsets[i]
        barplot!(barpositions, coverage;
            width=bar_width,
            label=est,
        )
    end

    axislegend(ax; position=:rt)
end

fig
