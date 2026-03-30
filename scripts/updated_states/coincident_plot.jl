using DrWatson
@quickactivate :PseudoStructuralComovements
using CSV, Tables, SparseArrays, Distributions, XLSX
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

state_names = ["IA", "IL", "IN", "MI", "MN", "ND", "OH", "SD", "WI"]
rawdata = XLSX.readdata(datadir("./state_indexes/reguib_northcentral.xlsx"), "Sheet1!A2:S459")
vecdata = Float64.(rawdata[:, 2:end])'
coincident = demean_standardize(vecdata[1:2:end, 278:end-2], dims = 2)

ps_cis = load(datadir("updated_states/coincident_series.jld2"))
ps_coincident = demean_standardize(ps_cis["cis"], dims = 2)
cor(ps_coincident[1,:], coincident[1,:])
cor(ps_coincident[2,:], coincident[2,:])
cor(ps_coincident[3,:], coincident[3,:])
cor(ps_coincident[4,:], coincident[4,:])
cor(ps_coincident[5,:], coincident[5,:])
cor(ps_coincident[6,:], coincident[6,:])
cor(ps_coincident[7,:], coincident[7,:])
cor(ps_coincident[8,:], coincident[8,:])
cor(ps_coincident[9,:], coincident[9,:])

state_names = ["IA", "IL", "IN", "MI", "MN", "ND", "OH", "SD", "WI"]
full_state_names = ["Iowa", "Illinois", "Indiana", "Michigan", "Minnesota", 
                    "North Dakota", "Ohio", "South Dakota", "Wisconsin"]

# Create a 3x3 grid of plots
plots = []
for i in 1:9
    p = Plots.plot(
        coincident[i, :], 
        label = "CCM",
        linewidth = 1.5,
        color = RGB(0/255, 114/255, 178/255),
        title = full_state_names[i],
        titlefontsize = 12,
        legend = i == 9 ? :bottomright : false,
        legendfontsize = 9,
        foreground_color_legend = nothing,  # removes the legend box border
        background_color_legend = RGBA(1,1,1,0.5),  # semi-transparent background
        xlabel = "",
        ylabel = "",
        xticks = (1:24:179, string.(2005:2:2019)),  # adjust to your time range
        xrotation = 45,
        tickfontsize = 6
    )
    Plots.plot!(p, 
        ps_coincident[i, :], 
        label = "Pseudo-structural",
        linewidth = 1.5,
        color = RGB(213/255, 94/255, 0/255),
    )
    push!(plots, p)
end

fig = Plots.plot(plots..., 
    layout = (3, 3), 
    size = (900, 700),
    margin = 3Plots.mm
)

savefig(fig, "coincident_indexes.pdf")
