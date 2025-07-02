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

smallbic31_bench = fill(NaN, 2, sims)
medbic31_bench = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin)
    med_bench_data = reshape(medmar.data', (medobs, dimvals[1], dimvals[2]))
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin)
    small_bench_data = reshape(smallmar.data', (smallobs, dimvals[1], dimvals[2]))

    smallicest = rank_selection(smallmar.data, dimvals; iters=100)
    smallaic31[:, s] .= smallicest.aic_sel[1:2]
    smallbic31[:, s] .= smallicest.bic_sel[1:2]

    medicest = rank_selection(medmar.data, dimvals; iters=100)
    medaic31[:, s] .= medicest.aic_sel[1:2]
    medbic31[:, s] .= medicest.bic_sel[1:2]

    small_bench = R"""
    d1 = $dimvals[1]
    d2 = $dimvals[2]
    small_data <- $small_bench_data
    small_selected_rank <- r_rank_selection(small_data, d1, d2)
    """
    @rget small_selected_rank

    med_bench = R"""
    d1 = $dimvals[1]
    d2 = $dimvals[2]
    med_data <- $med_bench_data
    med_selected_rank <- r_rank_selection(med_data, d1, d2)
    """
    @rget med_selected_rank

    smallbic31_bench[:, s] .= small_selected_rank[:selected_ranks]
    medbic31_bench[:, s] .= med_selected_rank[:selected_ranks]
end

save(datadir("threebyfour/31_results.jld2"), Dict(
    "smallaic" => smallaic31,
    "smallbic" => smallbic31,
    "medaic" => medaic31,
    "medbic" => medbic31,
    "smallbic_bench" => smallbic31_bench,
    "medbic_bench" => medbic31_bench,
))

medaicstats = sim_stats(medaic31, ranks, sims)
medbicstats = sim_stats(medbic31, ranks, sims)
smallaicstats = sim_stats(smallaic31, ranks, sims)
smallbicstats = sim_stats(smallbic31, ranks, sims)
smallbicstats_bench = sim_stats(smallbic31_bench, ranks, sims)
medbicstats_bench = sim_stats(medbic31_bench, ranks, sims)

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
