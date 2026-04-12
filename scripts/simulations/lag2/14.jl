using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250723)

sims = 100
dimvals = [3, 4]
ranks = [1, 4]
snr = 0.7

smallaic14 = fill(NaN, 3, sims)
smallbic14 = fill(NaN, 3, sims)
smallebic14 = fill(NaN, 3, sims)

medaic14 = fill(NaN, 3, sims)
medbic14 = fill(NaN, 3, sims)
medebic14 = fill(NaN, 3, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks; p=2)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin, matrix_err=true, p=2)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin, matrix_err=true, p=2)

    smallicest = rank_selection(smallmar.data, dimvals; iters=1000, pmax=2, num_starts=20, num_selected=3)
    smallaic14[:, s] .= smallicest.aic_sel[1:3]
    smallbic14[:, s] .= smallicest.bic_sel[1:3]
    smallebic14[:, s] .= smallicest.ebic_sel[1:3]

    medicest = rank_selection(medmar.data, dimvals; iters=1000, pmax=2, num_starts=20, num_selected=3)
    medaic14[:, s] .= medicest.aic_sel[1:3]
    medbic14[:, s] .= medicest.bic_sel[1:3]
    medebic14[:, s] .= medicest.ebic_sel[1:3]

end

save(datadir("lag2/14_results.jld2"), Dict(
    "smallaic" => smallaic14,
    "smallbic" => smallbic14,
    "smallebic" => smallebic14,
    "medaic" => medaic14,
    "medbic" => medbic14,
    "medebic" => medebic14,
))

push!(ranks, 1)
medaicstats = sim_stats(medaic14, ranks, sims)
medbicstats = sim_stats(medbic14, ranks, sims)
medebicstats = sim_stats(medebic14, ranks, sims)

smallaicstats = sim_stats(smallaic14, ranks, sims)
smallbicstats = sim_stats(smallbic14, ranks, sims)
smallebicstats = sim_stats(smallebic14, ranks, sims)

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
