using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250607)

sims = 100
dimvals = [8, 3]
ranks = [8, 1]
snr = 0.7

smallaic81 = fill(NaN, 2, sims)
smallbic81 = fill(NaN, 2, sims)
medaic81 = fill(NaN, 2, sims)
medbic81 = fill(NaN, 2, sims)

under_smallaic81 = fill(NaN, 2, sims)
under_smallbic81 = fill(NaN, 2, sims)
under_medaic81 = fill(NaN, 2, sims)
under_medbic81 = fill(NaN, 2, sims)

second_smallaic81 = fill(NaN, 2, sims)
second_smallbic81 = fill(NaN, 2, sims)
second_medaic81 = fill(NaN, 2, sims)
second_medbic81 = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin)

    smallicest = rank_selection(smallmar.data, dimvals; iters=200)
    smallaic81[:, s] .= smallicest.aic_sel[1:2]
    smallbic81[:, s] .= smallicest.bic_sel[1:2]

    small_clipped_icest = check_under_rank(smallicest.ictable, 2.0)
    under_smallaic81[:, s] .= small_clipped_icest.aic_sel[1:2]
    under_smallbic81[:, s] .= small_clipped_icest.bic_sel[1:2]

    second_small_clipped_icest = check_under_rank(smallicest.ictable, 1.0)
    second_smallaic81[:, s] .= second_small_clipped_icest.aic_sel[1:2]
    second_smallbic81[:, s] .= second_small_clipped_icest.bic_sel[1:2]

    ############################################################################

    medicest = rank_selection(medmar.data, dimvals; iters=200)
    medaic81[:, s] .= medicest.aic_sel[1:2]
    medbic81[:, s] .= medicest.bic_sel[1:2]

    med_clipped_icest = check_under_rank(medicest.ictable, 2.0)
    under_medaic81[:, s] .= med_clipped_icest.aic_sel[1:2]
    under_medbic81[:, s] .= med_clipped_icest.bic_sel[1:2]

    second_med_clipped_icest = check_under_rank(medicest.ictable, 1.0)
    second_medaic81[:, s] .= second_med_clipped_icest.aic_sel[1:2]
    second_medbic81[:, s] .= second_med_clipped_icest.bic_sel[1:2]
end

save(datadir("eightbythree/81_results.jld2"), Dict(
    "smallaic" => smallaic81,
    "smallbic" => smallbic81,
    "medaic" => medaic81,
    "medbic" => medbic81
))

save(datadir("eightbythree/under_81_results.jld2"), Dict(
    "under_smallaic" => under_smallaic81,
    "under_smallbic" => under_smallbic81,
    "under_medaic" => under_medaic81,
    "under_medbic" => under_medbic81
))

save(datadir("eightbythree/second_81_results.jld2"), Dict(
    "second_smallaic" => second_smallaic81,
    "second_smallbic" => second_smallbic81,
    "second_medaic" => second_medaic81,
    "second_medbic" => second_medbic81
))

medaicstats = sim_stats(medaic81, ranks, sims)
medbicstats = sim_stats(medbic81, ranks, sims)
smallaicstats = sim_stats(smallaic81, ranks, sims)
smallbicstats = sim_stats(smallbic81, ranks, sims)

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

under_medaicstats = sim_stats(under_medaic81, [2, 1], sims)
under_medbicstats = sim_stats(under_medbic81, [2, 1], sims)
under_smallaicstats = sim_stats(under_smallaic81, [2, 1], sims)
under_smallbicstats = sim_stats(under_smallbic81, [2, 1], sims)

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

second_medaicstats = sim_stats(second_medaic81, [1, 1], sims)
second_medbicstats = sim_stats(second_medbic81, [1, 1], sims)
second_smallaicstats = sim_stats(second_smallaic81, [1, 1], sims)
second_smallbicstats = sim_stats(second_smallbic81, [1, 1], sims)

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
