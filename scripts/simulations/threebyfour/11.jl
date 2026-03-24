using DrWatson
@quickactivate :PseudoStructuralComovements
using RCall
Random.seed!(20250723)

sims = 100
dimvals = [3, 4]
ranks = [1, 1]
snr = 0.7

smallaic11 = fill(NaN, 2, sims)
smallbic11 = fill(NaN, 2, sims)
smallebic11 = fill(NaN, 2, sims)
smallbench11 = fill(NaN, 2, sims)

medaic11 = fill(NaN, 2, sims)
medbic11 = fill(NaN, 2, sims)
medebic11 = fill(NaN, 2, sims)
medbench11 = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin, matrix_err=true)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin, matrix_err=true)

    smallicest = rank_selection(smallmar.data, dimvals; iters=1000)
    smallbench = rrmar_ic(smallmar.data, dimvals)
    smallaic11[:, s] .= smallicest.aic_sel[1:2]
    smallbic11[:, s] .= smallicest.bic_sel[1:2]
    smallebic11[:, s] .= smallicest.ebic_sel[1:2]
    smallbench11[:, s] .= smallbench.ebic_sel[1:2]

    medicest = rank_selection(medmar.data, dimvals; iters=1000)
    medbench = rrmar_ic(medmar.data, dimvals)
    medaic11[:, s] .= medicest.aic_sel[1:2]
    medbic11[:, s] .= medicest.bic_sel[1:2]
    medebic11[:, s] .= medicest.ebic_sel[1:2]
    medbench11[:, s] .= medbench.ebic_sel[1:2]

end

# save(datadir("threebyfour/11_results.jld2"), Dict(
#     "smallaic" => smallaic11,
#     "smallbic" => smallbic11,
#     "medaic" => medaic11,
#     "medbic" => medbic11,
#     "smallbic_bench" => smallbic11_bench,
#     "medbic_bench" => medbic11_bench,
# ))

medaicstats = sim_stats(medaic11, ranks, sims)
medbicstats = sim_stats(medbic11, ranks, sims)
smallaicstats = sim_stats(smallaic11, ranks, sims)
smallbicstats = sim_stats(smallbic11, ranks, sims)
smallbicstats_bench = sim_stats(smallbic11_bench, ranks, sims)
medbicstats_bench = sim_stats(medbic11_bench, ranks, sims)

println("Average rank for small size (AIC): ", smallaicstats.mad)
println("Average rank for small size (BIC): ", smallbicstats.mad)
println("Average rank for small size (bench BIC): ", smallbicstats_bench.mad)

println("Std. Dev rank for small size (AIC): ", round.(smallaicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(smallbicstats.stdval, digits=4))
println("Std. Dev rank for small size (bench BIC): ", round.(smallbicstats_bench.stdval, digits=4))

println("Freq. Correct for small size (AIC): ", smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)
println("Freq. Correct for small size (bench BIC): ", smallbicstats_bench.freqcorrect)

println("Average rank for medium size (AIC): ", medaicstats.mad)
println("Average rank for medium size (BIC): ", medbicstats.mad)
println("Average rank for medium size (bench BIC): ", medbicstats_bench.mad)

println("Std. Dev rank for medium size (AIC): ", round.(medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(medbicstats.stdval, digits=4))
println("Std. Dev rank for medium size (bench BIC): ", round.(medbicstats_bench.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
println("Freq. Correct for medium size (bench BIC): ", medbicstats_bench.freqcorrect)
