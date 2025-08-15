using DrWatson
@quickactivate :PseudoStructuralComovements
using Dates, XLSX

predata = XLSX.readdata(datadir("./state_indexes/reguib_northcentral.xlsx"), "Sheet1!A2:S459")
vecdata = Float64.(predata[:, 2:end])'
matdata = reshape(vecdata, (2,9,size(vecdata, 2)))
N1, N2, obs = size(matdata)

# build the date vector (for labels only)
start = Date(1982,1,1)
dates = start .+ Month(1)*(0:(obs-1))

# pick every 20th quarter (i.e. every 5 years) for ticks
tick_ix = 1:96:obs
tick_labels = lpad.(mod.(year.(dates[tick_ix]), 100), 2, "0")

states = ["IA", "IL", "IN", "MI", "MN", "ND", "OH", "SD", "WI"]
econ_indicators = ["CI", "LI"]

fig = Figure(backgroundcolor=:transparent, size=(1200, 350));
for i in 1:N1, j in 1:N2
    if (i == N1)
        ax = Axis(fig[i, j];
            xticks = (tick_ix, tick_labels),
            xticklabelrotation = π/6,
            xticklabelsize = 10
        )
    else
        ax = Axis(fig[i, j];
            xticks = (tick_ix, repeat([""], size(tick_labels, 1))),
            #=xticks = (tick_ix, tick_labels),=#
            #=xticklabelrotation = π/6=#
        )
    end

    lines!(ax, 1:obs, matdata[i, j, :])
    ax.title = "$(states[j]) $(econ_indicators[i])"
end

save("allstateplots.png", fig; transparent = true)

