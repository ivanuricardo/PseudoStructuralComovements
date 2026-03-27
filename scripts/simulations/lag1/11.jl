using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250723)

sims = 100
dimvals = [3, 4]
ranks = [1, 1]
snr = 0.7

smallaic11 = fill(NaN, 2, sims)
smallbic11 = fill(NaN, 2, sims)
smallebic11 = fill(NaN, 2, sims)

medaic11 = fill(NaN, 2, sims)
medbic11 = fill(NaN, 2, sims)
medebic11 = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin, matrix_err=true)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin, matrix_err=true)

    smallicest = rank_selection(smallmar.data, dimvals; iters=1000, pmax=2, num_starts=10, num_selected=4)
    smallaic11[:, s] .= smallicest.aic_sel[1:2]
    smallbic11[:, s] .= smallicest.bic_sel[1:2]
    smallebic11[:, s] .= smallicest.ebic_sel[1:2]

    medicest = rank_selection(medmar.data, dimvals; iters=1000, pmax=2, num_starts=10, num_selected=4)
    medaic11[:, s] .= medicest.aic_sel[1:2]
    medbic11[:, s] .= medicest.bic_sel[1:2]
    medebic11[:, s] .= medicest.ebic_sel[1:2]

end

save(datadir("lag1/11_results.jld2"), Dict(
    "smallaic" => smallaic11,
    "smallbic" => smallbic11,
    "smallebic" => smallebic11,
    "medaic" => medaic11,
    "medbic" => medbic11,
    "medebic" => medebic11,
))

push!(ranks, 1)
medaicstats = sim_stats(medaic11, ranks, sims)
medbicstats = sim_stats(medbic11, ranks, sims)
medebicstats = sim_stats(medebic11, ranks, sims)

smallaicstats = sim_stats(smallaic11, ranks, sims)
smallbicstats = sim_stats(smallbic11, ranks, sims)
smallebicstats = sim_stats(smallebic11, ranks, sims)

println("MAD for small size (AIC): ", smallaicstats.mad)
println("MAD for small size (BIC): ", smallbicstats.mad)
println("MAD for small size (EBIC): ", smallebicstats.mad)

println("Freq. Correct for small size (AIC): ", smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)
println("Freq. Correct for small size (EBIC): ", smallebicstats.freqcorrect)

println("MAD for medium size (AIC): ", medaicstats.mad)
println("MAD for medium size (BIC): ", medbicstats.mad)
println("MAD for medium size (EBIC): ", medebicstats.mad)

println("Freq. Correct for medium size (AIC): ", medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
println("Freq. Correct for medium size (EBIC): ", medebicstats.freqcorrect)
