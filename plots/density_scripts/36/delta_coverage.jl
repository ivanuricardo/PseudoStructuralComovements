using DrWatson
@quickactivate :PseudoStructuralComovements

using CairoMakie

colorblind_palette = [
    RGB(230/255, 159/255, 0/255),     # Orange
    RGB(86/255, 180/255, 233/255),    # Sky blue
    RGB(0/255, 158/255, 115/255),     # Bluish green
    RGB(240/255, 228/255, 66/255),    # Yellow
    RGB(0/255, 114/255, 178/255),     # Blue
    RGB(213/255, 94/255, 0/255),      # Vermillion
    RGB(204/255, 121/255, 167/255),   # Reddish purple
    RGB(0/255, 0/255, 0/255),         # Black
]

density_results100 = load(datadir("coverage/36/delta_cov_results100.jld2"))
density_results250 = load(datadir("coverage/36/delta_cov_results250.jld2"))

correct_cov100 = density_results100["correct_cov"]
under_cov100 = density_results100["under_cov"]
correct100 = vec(sum(correct_cov100, dims=2))
under100 = vec(sum(under_cov100, dims=2))

correct_cov250 = density_results250["correct_cov"]
under_cov250 = density_results250["under_cov"]
correct250 = vec(sum(correct_cov250, dims=2))
under250 = vec(sum(under_cov250, dims=2))

# Labels
deltas = ["δ₁", "δ₂"]
gamma_ranks = ["Correct", "Underestimated"]
sample_sizes = ["T = 100", "T = 250"]

# Now build coverage_data by looping over the three γ’s:
coverage_data = Dict{String,Vector{Vector{Float64}}}()
for (i, γ) in enumerate(deltas)
    coverage_data[γ] = [
        [correct100[i], under100[i]],  # T = 100
        [correct250[i], under250[i]]   # T = 250
    ]
end

# Convert to percentage
coverage_data = Dict(k => [map(x -> x / 10, row) for row in v] for (k, v) in coverage_data)

# Prepare plotting data
bar_width = 0.25
x_positions = [1, 2]  # Delta rank groups
offsets = [-bar_width / 2, bar_width / 2]  # For each gamma estimator

fig = Figure(
    size=(800, 500),
    backgroundcolor=:transparent     # make the figure background invisible
);

for (j, T) in enumerate(sample_sizes)
    ax = Axis(fig[j, 1];
        title=T,
        ylabel="Coverage (%)",
        xticks=(x_positions, gamma_ranks),
        limits=(nothing, (80, 100)),    # now from 80% up to 100%
        yticks=80:5:100,                # ticks every 5%
        backgroundcolor=:transparent,
    )

    for (i, est) in enumerate(deltas)
        coverage = coverage_data[est][j]
        barpositions = x_positions .+ offsets[i]
        barplot!(barpositions, coverage;
            width=bar_width,
            label=est,
            color=colorblind_palette[i],
        )
    end

    if j == 1
        axislegend(ax; position=:rt, backgroundcolor=:transparent)
    end
end

fig

save("../../density_plots/36/delta_coverage_plot_36.png", fig)
