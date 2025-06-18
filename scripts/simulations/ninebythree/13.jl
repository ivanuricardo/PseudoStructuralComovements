using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250607)

sims = 100
dimvals = [8, 3]
ranks = [1, 3]
snr = 0.7

smallaic13 = fill(NaN, 2, sims)
smallbic13 = fill(NaN, 2, sims)
medaic13 = fill(NaN, 2, sims)
medbic13 = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin)

    smallicest = rank_selection(smallmar.data, dimvals; iters=200)
    smallaic13[:, s] .= smallicest.aic_sel[1:2]
    smallbic13[:, s] .= smallicest.bic_sel[1:2]

    medicest = rank_selection(medmar.data, dimvals; iters=200)
    medaic13[:, s] .= medicest.aic_sel[1:2]
    medbic13[:, s] .= medicest.bic_sel[1:2]
end

save(datadir("eightbythree/13_results.jld2"), Dict(
    "smallaic" => smallaic13,
    "smallbic" => smallbic13,
    "medaic" => medaic13,
    "medbic" => medbic13
))

medaicstats = sim_stats(medaic13, ranks, sims)
medbicstats = sim_stats(medbic13, ranks, sims)
smallaicstats = sim_stats(smallaic13, ranks, sims)
smallbicstats = sim_stats(smallbic13, ranks, sims)

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
