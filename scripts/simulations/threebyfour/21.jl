using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250723)

sims = 100
dimvals = [3, 4]
ranks = [2, 1]
snr = 0.7

smallaic21 = fill(NaN, 2, sims)
smallbic21 = fill(NaN, 2, sims)
smallebic21 = fill(NaN, 2, sims)
smallbench21 = fill(NaN, 2, sims)

medaic21 = fill(NaN, 2, sims)
medbic21 = fill(NaN, 2, sims)
medebic21 = fill(NaN, 2, sims)
medbench21 = fill(NaN, 2, sims)

under_smallaic21 = fill(NaN, 2, sims)
under_smallbic21 = fill(NaN, 2, sims)
under_smallebic21 = fill(NaN, 2, sims)
under_smallbench21 = fill(NaN, 2, sims)

under_medaic21 = fill(NaN, 2, sims)
under_medbic21 = fill(NaN, 2, sims)
under_medebic21 = fill(NaN, 2, sims)
under_medbench21 = fill(NaN, 2, sims)

over_smallaic21 = fill(NaN, 2, sims)
over_smallbic21 = fill(NaN, 2, sims)
over_smallebic21 = fill(NaN, 2, sims)
over_smallbench21 = fill(NaN, 2, sims)

over_medaic21 = fill(NaN, 2, sims)
over_medbic21 = fill(NaN, 2, sims)
over_medebic21 = fill(NaN, 2, sims)
over_medbench21 = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin, matrix_err=true)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin, matrix_err=true)

    smallicest = rank_selection(smallmar.data, dimvals; iters=1000, num_starts=20, num_selected=3)
    smallbench = rrmar_ic(smallmar.data, dimvals)
    smallaic21[:, s] .= smallicest.aic_sel[1:2]
    smallbic21[:, s] .= smallicest.bic_sel[1:2]
    smallebic21[:, s] .= smallicest.ebic_sel[1:2]
    smallbench21[:, s] .= smallbench.ebic_sel[1:2]

    small_clipped_icest = check_rank(smallicest.ictable, 1.0)
    small_clipped_bench = check_rank(smallbench.ictable, 1.0)
    under_smallaic21[:, s] .= small_clipped_icest.aic_sel[1:2]
    under_smallbic21[:, s] .= small_clipped_icest.bic_sel[1:2]
    under_smallebic21[:, s] .= small_clipped_icest.ebic_sel[1:2]
    under_smallbench21[:, s] .= small_clipped_bench.ebic_sel[1:2]

    over_small_clipped_icest = check_rank(smallicest.ictable, 3.0)
    over_small_clipped_bench = check_rank(smallbench.ictable, 3.0)
    over_smallaic21[:, s] .= over_small_clipped_icest.aic_sel[1:2]
    over_smallbic21[:, s] .= over_small_clipped_icest.bic_sel[1:2]
    over_smallebic21[:, s] .= over_small_clipped_icest.ebic_sel[1:2]
    over_smallbench21[:, s] .= over_small_clipped_bench.ebic_sel[1:2]

    ############################################################################

    medicest = rank_selection(medmar.data, dimvals; iters=1000, num_starts=20, num_selected=3)
    medbench = rrmar_ic(medmar.data, dimvals)
    medaic21[:, s] .= medicest.aic_sel[1:2]
    medbic21[:, s] .= medicest.bic_sel[1:2]
    medebic21[:, s] .= medicest.ebic_sel[1:2]
    medbench21[:, s] .= medbench.ebic_sel[1:2]

    med_clipped_icest = check_rank(medicest.ictable, 1.0)
    med_clipped_bench = check_rank(medbench.ictable, 1.0)
    under_medaic21[:, s] .= med_clipped_icest.aic_sel[1:2]
    under_medbic21[:, s] .= med_clipped_icest.bic_sel[1:2]
    under_medebic21[:, s] .= med_clipped_icest.ebic_sel[1:2]
    under_medbench21[:, s] .= med_clipped_bench.ebic_sel[1:2]

    over_med_clipped_icest = check_rank(medicest.ictable, 3.0)
    over_med_clipped_bench = check_rank(medbench.ictable, 3.0)
    over_medaic21[:, s] .= over_med_clipped_icest.aic_sel[1:2]
    over_medbic21[:, s] .= over_med_clipped_icest.bic_sel[1:2]
    over_medebic21[:, s] .= over_med_clipped_icest.ebic_sel[1:2]
    over_medbench21[:, s] .= over_med_clipped_bench.ebic_sel[1:2]

end

save(datadir("threebyfour/21_results.jld2"), Dict(
    "smallaic" => smallaic21,
    "smallbic" => smallbic21,
    "smallebic" => smallebic21,
    "smallbench" => smallbench21,
    "medaic" => medaic21,
    "medbic" => medbic21,
    "medebic" => medebic21,
    "medbench" => medbench21,
))

save(datadir("threebyfour/under_21_results.jld2"), Dict(
    "under_smallaic" => under_smallaic21,
    "under_smallbic" => under_smallbic21,
    "under_smallebic" => under_smallebic21,
    "under_smallbench" => under_smallbench21,
    "under_medaic" => under_medaic21,
    "under_medbic" => under_medbic21,
    "under_medebic" => under_medebic21,
    "under_medbench" => under_medbench21,
))

save(datadir("threebyfour/over_21_results.jld2"), Dict(
    "over_smallaic" => over_smallaic21,
    "over_smallbic" => over_smallbic21,
    "over_smallebic" => over_smallebic21,
    "over_smallbench" => over_smallbench21,
    "over_medaic" => over_medaic21,
    "over_medbic" => over_medbic21,
    "over_medebic" => over_medebic21,
    "over_medbench" => over_medbench21,
))

medaicstats = sim_stats(medaic21, ranks, sims)
medbicstats = sim_stats(medbic21, ranks, sims)
medebicstats = sim_stats(medebic21, ranks, sims)
medbenchstats = sim_stats(medbench21, ranks, sims)

smallaicstats = sim_stats(smallaic21, ranks, sims)
smallbicstats = sim_stats(smallbic21, ranks, sims)
smallebicstats = sim_stats(smallebic21, ranks, sims)
smallbenchstats = sim_stats(smallbench21, ranks, sims)

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
println("Freq. Correct for medium size (Bench EBIC): ", medbenchstats.freqcorrect)

println("________________________________________________")

println("BONUS: Statistics if we 'accidentally' choose the incorrect rank")

under_medaicstats = sim_stats(under_medaic21, ranks, sims)
under_medbicstats = sim_stats(under_medbic21, ranks, sims)
under_medebicstats = sim_stats(under_medebic21, ranks, sims)
under_medbenchstats = sim_stats(under_medbench21, ranks, sims)

under_smallaicstats = sim_stats(under_smallaic21, ranks, sims)
under_smallbicstats = sim_stats(under_smallbic21, ranks, sims)
under_smallebicstats = sim_stats(under_smallebic21, ranks, sims)
under_smallbenchstats = sim_stats(under_smallbench21, ranks, sims)

println("MAD for small size (AIC): ", under_smallaicstats.mad)
println("MAD for small size (BIC): ", under_smallbicstats.mad)
println("MAD for small size (EBIC): ", under_smallebicstats.mad)
println("MAD for small size (Bench EBIC): ", under_smallbenchstats.mad)

println("Freq. Correct for small size (AIC): ", under_smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", under_smallbicstats.freqcorrect)
println("Freq. Correct for small size (EBIC): ", under_smallebicstats.freqcorrect)
println("Freq. Correct for small size (Bench EBIC): ", under_smallbenchstats.freqcorrect)

println("MAD for medium size (AIC): ", under_medaicstats.mad)
println("MAD for medium size (BIC): ", under_medbicstats.mad)
println("MAD for medium size (EBIC): ", under_medebicstats.mad)
println("MAD for medium size (Bench EBIC): ", under_medbenchstats.mad)

println("Freq. Correct for medium size (AIC): ", under_medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", under_medbicstats.freqcorrect)
println("Freq. Correct for medium size (EBIC): ", under_medebicstats.freqcorrect)
println("Freq. Correct for medium size (Bench EBIC): ", under_medbenchstats.freqcorrect)

println("________________________________________________")

over_medaicstats = sim_stats(over_medaic21, ranks, sims)
over_medbicstats = sim_stats(over_medbic21, ranks, sims)
over_medebicstats = sim_stats(over_medebic21, ranks, sims)
over_medbenchstats = sim_stats(over_medbench21, ranks, sims)

over_smallaicstats = sim_stats(over_smallaic21, ranks, sims)
over_smallbicstats = sim_stats(over_smallbic21, ranks, sims)
over_smallebicstats = sim_stats(over_smallebic21, ranks, sims)
over_smallbenchstats = sim_stats(over_smallbench21, ranks, sims)

println("MAD for small size (AIC): ", over_smallaicstats.mad)
println("MAD for small size (BIC): ", over_smallbicstats.mad)
println("MAD for small size (EBIC): ", over_smallebicstats.mad)
println("MAD for small size (Bench EBIC): ", over_smallbenchstats.mad)

println("Freq. Correct for small size (AIC): ", over_smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", over_smallbicstats.freqcorrect)
println("Freq. Correct for small size (EBIC): ", over_smallebicstats.freqcorrect)
println("Freq. Correct for small size (Bench EBIC): ", over_smallbenchstats.freqcorrect)

println("MAD for medium size (AIC): ", over_medaicstats.mad)
println("MAD for medium size (BIC): ", over_medbicstats.mad)
println("MAD for medium size (EBIC): ", over_medebicstats.mad)
println("MAD for medium size (Bench EBIC): ", over_medbenchstats.mad)

println("Freq. Correct for medium size (AIC): ", over_medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", over_medbicstats.freqcorrect)
println("Freq. Correct for medium size (EBIC): ", over_medebicstats.freqcorrect)
println("Freq. Correct for medium size (Bench EBIC): ", over_medbenchstats.freqcorrect)
