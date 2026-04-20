using DrWatson
@quickactivate :PseudoStructuralComovements
using Tables, Dates, CSV, XLSX
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

ut_employment = load_series("employment")
ut_unemployment = load_series("unemployment")
ut_hours = load_series("hours")

employment = transform(ut_employment)
unemployment = transform(ut_unemployment; type = "diff")
hours = transform(ut_hours)

file = XLSX.readxlsx(datadir("updated_states/wages.xlsx"))
sheet = file["Table"]
ut_wages = Float64.(XLSX.getdata(sheet, "C7:BJ15"))'
# Rearrange to be alphabetical and align with other series
# iowa, illinois, indiana, michigan, minnesota, north dakota, ohio, south dakota, wisconsin
rearranged_wages = ut_wages[:, [3,1,2,4,5,6,7,8,9]]
monthly_wages = quarterly_to_monthly(rearranged_wages)
wages = transform(monthly_wages)

catted_data = cat(employment', unemployment', hours', wages'; dims = 3)
tendata = permutedims(catted_data, (1,3,2))
# Corresponds to the perm WI, ND, OH, MN, MI, IA, SD, IL, IN
# perm_states = [9, 6, 7, 5, 4, 1, 8, 2, 3]
# Corresponds to the perm ND, IN, MI, MN, OH, SD, WI, IL, IA
perm_states = [6, 3, 4, 5, 7, 8, 9, 2, 1]
rearranged_tendata = tendata[perm_states, :, :]
matdata = vectorize(rearranged_tendata)

ut_tendata = cat(ut_employment', ut_unemployment', ut_hours', monthly_wages', dims = 3)
new_tendata = permutedims(ut_tendata, (1,3,2))
N1, N2, obs = size(new_tendata)

# build the date vector (for labels only)
start = Date(2005,2,1)
dates = start .+ Month(1)*(0:(obs-1))

# pick every 20th quarter (i.e. every 5 years) for ticks
tick_ix = 1:48:obs
tick_labels = string.(year.(dates[tick_ix]), " M", ((month.(dates[tick_ix]) .- 1) .÷ 3) .+ 1)

states = ["IA", "IL", "IN", "MI", "MN", "ND", "OH", "SD", "WI"]
econ_indicators = ["EMP", "UNR", "HOURS", "WAGES"]

fig = Figure(backgroundcolor=:transparent, size=(1200, 1600));
for i in 1:N1, j in 1:N2
    if (i == N1)
        ax = Axis(fig[i, j];
            xticks = (tick_ix, tick_labels),
            xticklabelrotation = π/6,
            titlesize=20,
            xticklabelsize = 16,
        )
    else
        ax = Axis(fig[i, j];
            xticks = (tick_ix, repeat([""], length(tick_ix))),
            titlesize=20,
            xticklabelsize = 16,
            #=xticks = (tick_ix, tick_labels),=#
            #=xticklabelrotation = π/6=#
        )
    end

    lines!(ax, 1:obs, new_tendata[i, j, :])
    ax.title = "$(states[i]) $(econ_indicators[j])"
end

save("states_econinds.png", fig; transparent = true)

