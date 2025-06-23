using DrWatson
@quickactivate :PseudoStructuralComovements
using RCall
R"""
source("r_helpers.R")
"""
Random.seed!(20250607)

sims = 100
dimvals = [3, 4]
ranks = [3, 1]
snr = 0.7

smallaic31 = fill(NaN, 2, sims)
smallbic31 = fill(NaN, 2, sims)
medaic31 = fill(NaN, 2, sims)
medbic31 = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin)
    med_data = reshape(medmar.data', (medobs, dimvals[1], dimvals[2]))
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin)
    small_data = reshape(smallmar.data', (smallobs, dimvals[1], dimvals[2]))

    small_res = R"""
    d1 = $dimvals[1]
    d2 = $dimvals[2]
    small_data <- $small_data
    selected_rank <- r_rank_selection(small_data, d1, d2)
    """
    med_res = R"""
    d1 = $dimvals[1]
    d2 = $dimvals[2]
    med_data <- $med_data
    selected_rank <- r_rank_selection(med_data, d1, d2)
    """
    small_res = rcopy(small_res)
    med_res = rcopy(med_res)

    smallbic31[:, s] .= small_res
    medbic31[:, s] .= med_res
end

save(datadir("threebyfour/31_bench_results.jld2"), Dict(
    "smallbic" => smallbic31,
    "medbic" => medbic31
))

medbicstats = sim_stats(medbic31, ranks, sims)
smallbicstats = sim_stats(smallbic31, ranks, sims)

println("Average rank for small size (BIC): ", smallbicstats.avgval)
println("Std. Dev rank for small size (BIC): ", round.(smallbicstats.stdval, digits=4))
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)
println("Average rank for medium size (BIC): ", medbicstats.avgval)
println("Std. Dev rank for medium size (BIC): ", round.(medbicstats.stdval, digits=4))
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
