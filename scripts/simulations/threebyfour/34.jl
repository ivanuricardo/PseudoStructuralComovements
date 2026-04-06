using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250723)

sims = 100
dimvals = [3, 4]
ranks = [3, 4]
snr = 0.7

smallaic34 = fill(NaN, 2, sims)
smallbic34 = fill(NaN, 2, sims)
smallebic34 = fill(NaN, 2, sims)
smallbench34 = fill(NaN, 2, sims)

medaic34 = fill(NaN, 2, sims)
medbic34 = fill(NaN, 2, sims)
medebic34 = fill(NaN, 2, sims)
medbench34 = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin, matrix_err=true)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin, matrix_err=true)

    smallicest = rank_selection(smallmar.data, dimvals; iters=1000)
    smallbench = rrmar_ic(smallmar.data, dimvals)
    smallaic34[:, s] .= smallicest.aic_sel[1:2]
    smallbic34[:, s] .= smallicest.bic_sel[1:2]
    smallebic34[:, s] .= smallicest.ebic_sel[1:2]
    smallbench34[:, s] .= smallbench.ebic_sel[1:2]

    medicest = rank_selection(medmar.data, dimvals; iters=1000)
    medbench = rrmar_ic(medmar.data, dimvals)
    medaic34[:, s] .= medicest.aic_sel[1:2]
    medbic34[:, s] .= medicest.bic_sel[1:2]
    medebic34[:, s] .= medicest.ebic_sel[1:2]
    medbench34[:, s] .= medbench.ebic_sel[1:2]

end

save(datadir("threebyfour/34_results.jld2"), Dict(
    "smallaic" => smallaic34,
    "smallbic" => smallbic34,
    "smallebic" => smallebic34,
    "smallbench" => smallbench34,
    "medaic" => medaic34,
    "medbic" => medbic34,
    "medebic" => medebic34,
    "medbench" => medbench34,
))

medaicstats = sim_stats(medaic34, ranks, sims)
medbicstats = sim_stats(medbic34, ranks, sims)
medebicstats = sim_stats(medebic34, ranks, sims)
medbenchstats = sim_stats(medbench34, ranks, sims)

smallaicstats = sim_stats(smallaic34, ranks, sims)
smallbicstats = sim_stats(smallbic34, ranks, sims)
smallebicstats = sim_stats(smallebic34, ranks, sims)
smallbenchstats = sim_stats(smallbench34, ranks, sims)

println("MAD for small size (AIC): ", smallaicstats.mad)
println("MAD for small size (BIC): ", smallbicstats.mad)
println("MAD for small size (EBIC): ", smallebicstats.mad)
println("MAD for small size (Bench EBIC): ", smallbenchstats.mad)

println("Freq. Correct for small size (AIC): ", smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)
println("Freq. Correct for small size (EBIC): ", smallebicstats.freqcorrect)
println("Freq. Correct for small size (Bench EBIC): ", smallbenchstats.freqcorrect)

println("MAD for medium size (AIC): ", medaicstats.mad)
println("MAD for medium size (BIC): ", medbicstats.mad)
println("MAD for medium size (EBIC): ", medebicstats.mad)
println("MAD for medium size (Bench EBIC): ", medbenchstats.mad)

println("Freq. Correct for medium size (AIC): ", medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
println("Freq. Correct for medium size (EBIC): ", medebicstats.freqcorrect)
println("Freq. Correct for medium size (Bench BIC): ", medbenchstats.freqcorrect)
