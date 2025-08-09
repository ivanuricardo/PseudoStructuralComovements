using DrWatson
@quickactivate :PseudoStructuralComovements
using RCall
R"""
source("r_helpers.R")
"""
Random.seed!(20250723)

sims = 100
dimvals = [3, 4]
ranks = [3, 1]
snr = 0.7

smallaic31 = fill(NaN, 3, sims)
smallbic31 = fill(NaN, 3, sims)
medaic31 = fill(NaN, 3, sims)
medbic31 = fill(NaN, 3, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin, matrix_err=true)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin, matrix_err=true)

    smallicest = rank_selection(smallmar.data, dimvals; iters=1000, pmax=2)
    smallaic31[:, s] .= smallicest.aic_sel[1:3]
    smallbic31[:, s] .= smallicest.bic_sel[1:3]

    medicest = rank_selection(medmar.data, dimvals; iters=1000, pmax=2)
    medaic31[:, s] .= medicest.aic_sel[1:3]
    medbic31[:, s] .= medicest.bic_sel[1:3]
end

save(datadir("lag1/31_results.jld2"), Dict(
    "smallaic" => smallaic31,
    "smallbic" => smallbic31,
    "medaic" => medaic31,
    "medbic" => medbic31,
    "smallbic_bench" => smallbic31_bench,
    "medbic_bench" => medbic31_bench,
))

push!(ranks, 1)
medaicstats = sim_stats(medaic31, ranks, sims)
medbicstats = sim_stats(medbic31, ranks, sims)
smallaicstats = sim_stats(smallaic31, ranks, sims)
smallbicstats = sim_stats(smallbic31, ranks, sims)
smallbicstats_bench = sim_stats(smallbic31_bench, ranks, sims)
medbicstats_bench = sim_stats(medbic31_bench, ranks, sims)

println("Average rank for small size (AIC): ", smallaicstats.avgval)
println("Average lag for small size (AIC): ", smallaicstats.avgval)
println("Average rank for small size (BIC): ", smallbicstats.avgval)
println("Average lag for small size (BIC): ", smallbicstats.avgval)

println("Std. Dev rank for small size (AIC): ", round.(smallaicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(smallbicstats.stdval, digits=4))

println("Freq. Correct for small size (AIC): ", smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)

println("Average rank for medium size (AIC): ", medaicstats.avgval)
println("Average rank for medium size (BIC): ", medbicstats.avgval)

println("Std. Dev rank for medium size (AIC): ", round.(medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(medbicstats.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
