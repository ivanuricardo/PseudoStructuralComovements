using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250324)

sims = 100
dimvals = [3, 4]
ranks = [1, 1]
snr = 0.7

smallaic11 = fill(NaN, 2, sims)
smallbic11 = fill(NaN, 2, sims)
medaic11 = fill(NaN, 2, sims)
medbic11 = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin)
    small_data = medmar.data[:, 1:100]

    smallicest = rank_selection(small_data, dimvals; iters=100)
    smallaic11[:, s] .= smallicest.aic_sel[1:2]
    smallbic11[:, s] .= smallicest.bic_sel[1:2]

    medicest = rank_selection(medmar.data, dimvals; iters=100)
    medaic11[:, s] .= medicest.aic_sel[1:2]
    medbic11[:, s] .= medicest.bic_sel[1:2]
end

save(datadir("fourbyfive/11_results.jld2"), Dict(
    "smallaic" => smallaic11,
    "smallbic" => smallbic11,
    "medaic" => medaic11,
    "medbic" => medbic11
))

medaicstats = sim_stats(medaic11, ranks, sims)
medbicstats = sim_stats(medbic11, ranks, sims)
smallaicstats = sim_stats(smallaic11, ranks, sims)
smallbicstats = sim_stats(smallbic11, ranks, sims)

println("Average rank for small size (AIC): ", smallaicstats.avgval)
println("Average rank for small size (BIC): ", smallbicstats.avgval)

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
