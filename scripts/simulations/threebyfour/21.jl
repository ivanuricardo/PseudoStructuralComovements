using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250607)

sims = 100
dimvals = [3, 4]
ranks = [2, 1]
snr = 0.7

smallaic21 = fill(NaN, 2, sims)
smallbic21 = fill(NaN, 2, sims)
medaic21 = fill(NaN, 2, sims)
medbic21 = fill(NaN, 2, sims)

under_smallaic21 = fill(NaN, 2, sims)
under_smallbic21 = fill(NaN, 2, sims)
under_medaic21 = fill(NaN, 2, sims)
under_medbic21 = fill(NaN, 2, sims)

over_smallaic21 = fill(NaN, 2, sims)
over_smallbic21 = fill(NaN, 2, sims)
over_medaic21 = fill(NaN, 2, sims)
over_medbic21 = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin)

    smallicest = rank_selection(smallmar.data, dimvals; iters=200)
    smallaic21[:, s] .= smallicest.aic_sel[1:2]
    smallbic21[:, s] .= smallicest.bic_sel[1:2]

    small_clipped_icest = check_rank(smallicest.ictable, 1.0)
    under_smallaic21[:, s] .= small_clipped_icest.aic_sel[1:2]
    under_smallbic21[:, s] .= small_clipped_icest.bic_sel[1:2]

    over_small_clipped_icest = check_rank(smallicest.ictable, 3.0)
    over_smallaic21[:, s] .= over_small_clipped_icest.aic_sel[1:2]
    over_smallbic21[:, s] .= over_small_clipped_icest.bic_sel[1:2]

    ############################################################################

    medicest = rank_selection(medmar.data, dimvals; iters=200)
    medaic21[:, s] .= medicest.aic_sel[1:2]
    medbic21[:, s] .= medicest.bic_sel[1:2]

    med_clipped_icest = check_rank(medicest.ictable, 1.0)
    under_medaic21[:, s] .= med_clipped_icest.aic_sel[1:2]
    under_medbic21[:, s] .= med_clipped_icest.bic_sel[1:2]

    over_med_clipped_icest = check_rank(medicest.ictable, 3.0)
    over_medaic21[:, s] .= over_med_clipped_icest.aic_sel[1:2]
    over_medbic21[:, s] .= over_med_clipped_icest.bic_sel[1:2]
end

save(datadir("threebyfour/21_results.jld2"), Dict(
    "smallaic" => smallaic21,
    "smallbic" => smallbic21,
    "medaic" => medaic21,
    "medbic" => medbic21
))

save(datadir("threebyfour/under_21_results.jld2"), Dict(
    "under_smallaic" => under_smallaic21,
    "under_smallbic" => under_smallbic21,
    "under_medaic" => under_medaic21,
    "under_medbic" => under_medbic21
))

save(datadir("threebyfour/over_21_results.jld2"), Dict(
    "over_smallaic" => over_smallaic21,
    "over_smallbic" => over_smallbic21,
    "over_medaic" => over_medaic21,
    "over_medbic" => over_medbic21
))

medaicstats = sim_stats(medaic21, ranks, sims)
medbicstats = sim_stats(medbic21, ranks, sims)
smallaicstats = sim_stats(smallaic21, ranks, sims)
smallbicstats = sim_stats(smallbic21, ranks, sims)

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

println("________________________________________________")

println("BONUS: Statistics if we 'accidentally' choose the incorrect rank")

under_medaicstats = sim_stats(under_medaic21, ranks, sims)
under_medbicstats = sim_stats(under_medbic21, ranks, sims)
under_smallaicstats = sim_stats(under_smallaic21, ranks, sims)
under_smallbicstats = sim_stats(under_smallbic21, ranks, sims)

println("Average rank for small size (AIC): ", under_smallaicstats.avgval)
println("Average rank for small size (BIC): ", under_smallbicstats.avgval)

println("Std. Dev rank for small size (AIC): ", round.(under_smallaicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(under_smallbicstats.stdval, digits=4))

println("Freq. Correct for small size (AIC): ", under_smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", under_smallbicstats.freqcorrect)

println("Average rank for medium size (AIC): ", under_medaicstats.avgval)
println("Average rank for medium size (BIC): ", under_medbicstats.avgval)

println("Std. Dev rank for medium size (AIC): ", round.(under_medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(under_medbicstats.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", under_medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", under_medbicstats.freqcorrect)

println("________________________________________________")

over_medaicstats = sim_stats(over_medaic21, ranks, sims)
over_medbicstats = sim_stats(over_medbic21, ranks, sims)
over_smallaicstats = sim_stats(over_smallaic21, ranks, sims)
over_smallbicstats = sim_stats(over_smallbic21, ranks, sims)

println("Average rank for small size (AIC): ", over_smallaicstats.avgval)
println("Average rank for small size (BIC): ", over_smallbicstats.avgval)

println("Std. Dev rank for small size (AIC): ", round.(over_smallaicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(over_smallbicstats.stdval, digits=4))

println("Freq. Correct for small size (AIC): ", over_smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", over_smallbicstats.freqcorrect)

println("Average rank for medium size (AIC): ", over_medaicstats.avgval)
println("Average rank for medium size (BIC): ", over_medbicstats.avgval)

println("Std. Dev rank for medium size (AIC): ", round.(over_medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(over_medbicstats.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", over_medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", over_medbicstats.freqcorrect)
