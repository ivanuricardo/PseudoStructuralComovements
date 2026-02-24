using DrWatson
@quickactivate :PseudoStructuralComovements
using Downloads, CSV
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

url = "https://www2.census.gov/econ/bps/State/st0001c.txt"
Downloads.download(url, datadir("./updated_states/raw-permits/st0001c.txt")

