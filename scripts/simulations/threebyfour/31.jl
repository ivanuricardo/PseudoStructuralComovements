using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250607)

sims = 100
dimvals = [3, 4]
ranks = [3, 1]
snr = 0.7

smallaic11 = fill(NaN, 2, sims)
smallbic11 = fill(NaN, 2, sims)
medaic11 = fill(NaN, 2, sims)
medbic11 = fill(NaN, 2, sims)

under_smallaic11 = fill(NaN, 2, sims)
under_smallbic11 = fill(NaN, 2, sims)
under_medaic11 = fill(NaN, 2, sims)
under_medbic11 = fill(NaN, 2, sims)

second_smallaic11 = fill(NaN, 2, sims)
second_smallbic11 = fill(NaN, 2, sims)
second_medaic11 = fill(NaN, 2, sims)
second_medbic11 = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin)

    smallicest = rank_selection(smallmar.data, dimvals; iters=100)
    smallaic11[:, s] .= smallicest.aic_sel[1:2]
    smallbic11[:, s] .= smallicest.bic_sel[1:2]

    small_clipped_icest = check_under_rank(smallicest.ictable, 2.0)
    under_smallaic11[:, s] .= small_clipped_icest.aic_sel[1:2]
    under_smallbic11[:, s] .= small_clipped_icest.bic_sel[1:2]

    second_small_clipped_icest = check_under_rank(smallicest.ictable, 1.0)
    second_smallaic11[:, s] .= small_clipped_icest.aic_sel[1:2]
    second_smallbic11[:, s] .= small_clipped_icest.bic_sel[1:2]

    ############################################################################

    medicest = rank_selection(medmar.data, dimvals; iters=100)
    medaic11[:, s] .= medicest.aic_sel[1:2]
    medbic11[:, s] .= medicest.bic_sel[1:2]

    med_clipped_icest = check_under_rank(medicest.ictable, 2.0)
    under_medaic11[:, s] .= med_clipped_icest.aic_sel[1:2]
    under_medbic11[:, s] .= med_clipped_icest.bic_sel[1:2]

    second_med_clipped_icest = check_under_rank(medicest.ictable, 1.0)
    second_medaic11[:, s] .= med_clipped_icest.aic_sel[1:2]
    second_medbic11[:, s] .= med_clipped_icest.bic_sel[1:2]
end

save(datadir("threebyfour/31_results.jld2"), Dict(
    "smallaic" => smallaic11,
    "smallbic" => smallbic11,
    "medaic" => medaic11,
    "medbic" => medbic11
))

save(datadir("threebyfour/under_31_results.jld2"), Dict(
    "under_smallaic" => under_smallaic11,
    "under_smallbic" => under_smallbic11,
    "under_medaic" => under_medaic11,
    "under_medbic" => under_medbic11
))

save(datadir("threebyfour/second_31_results.jld2"), Dict(
    "second_smallaic" => second_smallaic11,
    "second_smallbic" => second_smallbic11,
    "second_medaic" => second_medaic11,
    "second_medbic" => second_medbic11
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

println("________________________________________________")

println("BONUS: Statistics if we 'accidentally' choose the incorrect rank")

under_medaicstats = sim_stats(under_medaic11, ranks, sims)
under_medbicstats = sim_stats(under_medbic11, ranks, sims)
under_smallaicstats = sim_stats(under_smallaic11, ranks, sims)
under_smallbicstats = sim_stats(under_smallbic11, ranks, sims)

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

second_medaicstats = sim_stats(second_medaic11, ranks, sims)
second_medbicstats = sim_stats(second_medbic11, ranks, sims)
second_smallaicstats = sim_stats(second_smallaic11, ranks, sims)
second_smallbicstats = sim_stats(second_smallbic11, ranks, sims)

println("Average rank for small size (AIC): ", second_smallaicstats.avgval)
println("Average rank for small size (BIC): ", second_smallbicstats.avgval)

println("Std. Dev rank for small size (AIC): ", round.(second_smallaicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(second_smallbicstats.stdval, digits=4))

println("Freq. Correct for small size (AIC): ", second_smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", second_smallbicstats.freqcorrect)

println("Average rank for medium size (AIC): ", second_medaicstats.avgval)
println("Average rank for medium size (BIC): ", second_medbicstats.avgval)

println("Std. Dev rank for medium size (AIC): ", round.(second_medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(second_medbicstats.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", second_medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", second_medbicstats.freqcorrect)
