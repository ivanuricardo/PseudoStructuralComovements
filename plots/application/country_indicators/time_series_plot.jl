using DrWatson
@quickactivate :PseudoStructuralComovements
using Dates

matdata = load(datadir("country_indicators/globaldata.jld2"), "matdata")
N1, N2, obs = size(matdata)

# build the date vector (for labels only)
start = Date(1991,1,1)
dates = start .+ Month(3)*(0:(obs-1))

# pick every 20th quarter (i.e. every 5 years) for ticks
tick_ix = 1:24:obs
tick_labels = string.(year.(dates[tick_ix]), " Q", ((month.(dates[tick_ix]) .- 1) .÷ 3) .+ 1)

countries = ["USA", "CAN", "DEU", "FRA"]
econ_indicators = ["GDP", "PROD", "IR"]

fig = Figure(backgroundcolor=:transparent, size=(800, 500));
for i in 1:N1, j in 1:N2
    if (i == N1)
        ax = Axis(fig[i, j];
            xticks = (tick_ix, tick_labels),
            xticklabelrotation = π/6
        )
    else
        ax = Axis(fig[i, j];
            xticks = (tick_ix, repeat([""], 5)),
            #=xticks = (tick_ix, tick_labels),=#
            #=xticklabelrotation = π/6=#
        )
    end

    lines!(ax, 1:obs, matdata[i, j, :])
    ax.title = "$(countries[j]) $(econ_indicators[i])"
end

save("country_indicators.png", fig; transparent = true)

