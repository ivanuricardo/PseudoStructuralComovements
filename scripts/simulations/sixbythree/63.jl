using DrWatson
@quickactivate :PseudoStructuralComovements
using RCall
R"""
source("r_helpers.R")
"""
Random.seed!(20250607)

sims = 100
dimvals = [6, 3]
ranks = [6, 3]
snr = 0.7

smallaic63 = fill(NaN, 2, sims)
smallbic63 = fill(NaN, 2, sims)
medaic63 = fill(NaN, 2, sims)
medbic63 = fill(NaN, 2, sims)

smallbic63_bench = fill(NaN, 2, sims)
medbic63_bench = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin)
    med_bench_data = reshape(medmar.data', (medobs, dimvals[1], dimvals[2]))
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin)
    small_bench_data = reshape(smallmar.data', (smallobs, dimvals[1], dimvals[2]))

    smallicest = rank_selection(smallmar.data, dimvals; iters=200)
    smallaic63[:, s] .= smallicest.aic_sel[1:2]
    smallbic63[:, s] .= smallicest.bic_sel[1:2]

    medicest = rank_selection(medmar.data, dimvals; iters=200)
    medaic63[:, s] .= medicest.aic_sel[1:2]
    medbic63[:, s] .= medicest.bic_sel[1:2]

    small_bench = R"""
    d1 = $dimvals[1]
    d2 = $dimvals[2]
    small_data <- $small_bench_data
    selected_rank <- r_rank_selection(small_data, d1, d2)
    """
    med_bench = R"""
    d1 = $dimvals[1]
    d2 = $dimvals[2]
    med_data <- $med_bench_data
    selected_rank <- r_rank_selection(med_data, d1, d2)
    """

    smallbic63_bench[:, s] .= rcopy(small_bench)
    medbic63_bench[:, s] .= rcopy(med_bench)
end

save(datadir("sixbythree/63_results.jld2"), Dict(
    "smallaic" => smallaic63,
    "smallbic" => smallbic63,
    "medaic" => medaic63,
    "medbic" => medbic63,
    "smallbic_bench" => smallbic63_bench,
    "medbic_bench" => medbic63_bench,
))

medaicstats = sim_stats(medaic63, ranks, sims)
medbicstats = sim_stats(medbic63, ranks, sims)
smallaicstats = sim_stats(smallaic63, ranks, sims)
smallbicstats = sim_stats(smallbic63, ranks, sims)
smallbicstats_bench = sim_stats(smallbic63_bench, ranks, sims)
medbicstats_bench = sim_stats(medbic63_bench, ranks, sims)

println("Average rank for small size (AIC): ", smallaicstats.avgval)
println("Average rank for small size (BIC): ", smallbicstats.avgval)
println("Average rank for small size (bench BIC): ", smallbicstats_bench.avgval)

println("Std. Dev rank for small size (AIC): ", round.(smallaicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(smallbicstats.stdval, digits=4))
println("Std. Dev rank for small size (bench BIC): ", round.(smallbicstats_bench.stdval, digits=4))

println("Freq. Correct for small size (AIC): ", smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)
println("Freq. Correct for small size (bench BIC): ", smallbicstats_bench.freqcorrect)

println("Average rank for medium size (AIC): ", medaicstats.avgval)
println("Average rank for medium size (BIC): ", medbicstats.avgval)
println("Average rank for medium size (bench BIC): ", medbicstats_bench.avgval)

println("Std. Dev rank for medium size (AIC): ", round.(medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(medbicstats.stdval, digits=4))
println("Std. Dev rank for medium size (bench BIC): ", round.(medbicstats_bench.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
println("Freq. Correct for medium size (bench BIC): ", medbicstats_bench.freqcorrect)
