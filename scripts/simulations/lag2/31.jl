using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250723)

sims = 100
dimvals = [3, 4]
ranks = [3, 1]
snr = 0.7

smallaic31 = fill(NaN, 3, sims)
smallbic31 = fill(NaN, 3, sims)
smallebic31 = fill(NaN, 3, sims)

medaic31 = fill(NaN, 3, sims)
medbic31 = fill(NaN, 3, sims)
medebic31 = fill(NaN, 3, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks; p=2)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin, matrix_err=true, p=2)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin, matrix_err=true, p=2)

    smallicest = rank_selection(smallmar.data, dimvals; iters=1000, pmax=2, num_starts=20, num_selected=3)
    smallaic31[:, s] .= smallicest.aic_sel[1:3]
    smallbic31[:, s] .= smallicest.bic_sel[1:3]
    smallebic31[:, s] .= smallicest.ebic_sel[1:3]

    medicest = rank_selection(medmar.data, dimvals; iters=1000, pmax=2, num_starts=20, num_selected=3)
    medaic31[:, s] .= medicest.aic_sel[1:3]
    medbic31[:, s] .= medicest.bic_sel[1:3]
    medebic31[:, s] .= medicest.ebic_sel[1:3]

end

save(datadir("lag2/31_results.jld2"), Dict(
    "smallaic" => smallaic31,
    "smallbic" => smallbic31,
    "smallebic" => smallebic31,
    "medaic" => medaic31,
    "medbic" => medbic31,
    "medebic" => medebic31,
))

push!(ranks, 1)
medaicstats = sim_stats(medaic31, ranks, sims)
medbicstats = sim_stats(medbic31, ranks, sims)
medebicstats = sim_stats(medebic31, ranks, sims)

smallaicstats = sim_stats(smallaic31, ranks, sims)
smallbicstats = sim_stats(smallbic31, ranks, sims)
smallebicstats = sim_stats(smallebic31, ranks, sims)

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
