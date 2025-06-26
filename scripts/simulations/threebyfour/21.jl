using DrWatson
@quickactivate :PseudoStructuralComovements
using RCall
R"""
source("r_helpers.R")
"""
Random.seed!(20250607)

sims = 3
dimvals = [3, 4]
ranks = [2, 1]
snr = 0.7

smallaic21 = fill(NaN, 2, sims)
smallbic21 = fill(NaN, 2, sims)
medaic21 = fill(NaN, 2, sims)
medbic21 = fill(NaN, 2, sims)

smallbic21_bench = fill(NaN, 2, sims)
medbic21_bench = fill(NaN, 2, sims)

under_smallaic21 = fill(NaN, 2, sims)
under_smallbic21 = fill(NaN, 2, sims)
under_medaic21 = fill(NaN, 2, sims)
under_medbic21 = fill(NaN, 2, sims)

under_smallbic21_bench = fill(NaN, 2, sims)
under_medbic21_bench = fill(NaN, 2, sims)

over_smallaic21 = fill(NaN, 2, sims)
over_smallbic21 = fill(NaN, 2, sims)
over_medaic21 = fill(NaN, 2, sims)
over_medbic21 = fill(NaN, 2, sims)

over_smallbic21_bench = fill(NaN, 2, sims)
over_medbic21_bench = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin)
    med_bench_data = reshape(medmar.data', (medobs, dimvals[1], dimvals[2]))
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin)
    small_bench_data = reshape(smallmar.data', (smallobs, dimvals[1], dimvals[2]))

    smallicest = rank_selection(smallmar.data, dimvals; iters=2)
    smallaic21[:, s] .= smallicest.aic_sel[1:2]
    smallbic21[:, s] .= smallicest.bic_sel[1:2]

    small_bench = R"""
    d1 = $dimvals[1]
    d2 = $dimvals[2]
    small_data <- $small_bench_data
    small_selected_rank <- r_rank_selection(small_data, d1, d2)
    """
    @rget small_selected_rank

    smallbic21_bench[:, s] .= small_selected_rank[:selected_ranks]

    small_clipped_icest = check_rank(smallicest.ictable, 1.0)
    under_smallaic21[:, s] .= small_clipped_icest.aic_sel[1:2]
    under_smallbic21[:, s] .= small_clipped_icest.bic_sel[1:2]

    under_small_bench_table = check_rank(small_selected_rank[:ic_results]', 1.0)
    under_smallbic21_bench[:, s] = under_small_bench_table.aic_sel

    over_small_clipped_icest = check_rank(smallicest.ictable, 3.0)
    over_smallaic21[:, s] .= over_small_clipped_icest.aic_sel[1:2]
    over_smallbic21[:, s] .= over_small_clipped_icest.bic_sel[1:2]

    over_small_bench_table = check_rank(small_selected_rank[:ic_results]', 3.0)
    over_smallbic21_bench[:, s] = over_small_bench_table.aic_sel

    ############################################################################

    medicest = rank_selection(medmar.data, dimvals; iters=2)
    medaic21[:, s] .= medicest.aic_sel[1:2]
    medbic21[:, s] .= medicest.bic_sel[1:2]

    med_bench = R"""
    d1 = $dimvals[1]
    d2 = $dimvals[2]
    med_data <- $med_bench_data
    med_selected_rank <- r_rank_selection(med_data, d1, d2)
    """
    @rget med_selected_rank

    medbic21_bench[:, s] .= med_selected_rank[:selected_ranks]

    med_clipped_icest = check_rank(medicest.ictable, 1.0)
    under_medaic21[:, s] .= med_clipped_icest.aic_sel[1:2]
    under_medbic21[:, s] .= med_clipped_icest.bic_sel[1:2]

    under_med_bench_table = check_rank(med_selected_rank[:ic_results]', 1.0)
    under_medbic21_bench[:, s] = under_med_bench_table.aic_sel

    over_med_clipped_icest = check_rank(medicest.ictable, 3.0)
    over_medaic21[:, s] .= over_med_clipped_icest.aic_sel[1:2]
    over_medbic21[:, s] .= over_med_clipped_icest.bic_sel[1:2]

    over_small_bench_table = check_rank(small_selected_rank[:ic_results]', 3.0)
    over_smallbic21_bench[:, s] = over_small_bench_table.aic_sel
end

save(datadir("threebyfour/21_results.jld2"), Dict(
    "smallaic" => smallaic21,
    "smallbic" => smallbic21,
    "smallbic_bench" => smallbic21_bench,
    "medaic" => medaic21,
    "medbic" => medbic21,
    "medbic_bench" => medbic21_bench,
))

save(datadir("threebyfour/under_21_results.jld2"), Dict(
    "under_smallaic" => under_smallaic21,
    "under_smallbic" => under_smallbic21,
    "under_smallbic_bench" => under_smallbic21_bench,
    "under_medaic" => under_medaic21,
    "under_medbic" => under_medbic21,
    "under_medbic_bench" => under_medbic21_bench,
))

save(datadir("threebyfour/over_21_results.jld2"), Dict(
    "over_smallaic" => over_smallaic21,
    "over_smallbic" => over_smallbic21,
    "over_smallbic_bench" => over_smallbic21_bench,
    "over_medaic" => over_medaic21,
    "over_medbic" => over_medbic21,
    "over_medbic_bench" => over_medbic21_bench,
))

medaicstats = sim_stats(medaic21, ranks, sims)
medbicstats = sim_stats(medbic21, ranks, sims)
medbicstats_bench = sim_stats(medbic21_bench, ranks, sims)
smallaicstats = sim_stats(smallaic21, ranks, sims)
smallbicstats = sim_stats(smallbic21, ranks, sims)
smallbicstats_bench = sim_stats(smallbic21_bench, ranks, sims)

println("Average rank for small size (AIC): ", smallaicstats.avgval)
println("Average rank for small size (BIC): ", smallbicstats.avgval)
println("Average rank for small size (bench BIC): ", smallbicstats_bench.avgval)

println("Std. Dev rank for small size (AIC): ", round.(smallaicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(smallbicstats.stdval, digits=4))
println("Std. Dev rank for small size (bench BIC): ", round.(smallbicstats_bench.stdval, digits=4))

println("Freq. Correct for small size (AIC): ", smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)
println("Freq. Correct for small size (bench BIC): ", smallbicstats_bench.freqcorrect)

println("Average rank for medium size (AIC): ", medaicstats.avgval)
println("Average rank for medium size (BIC): ", medbicstats.avgval)
println("Average rank for medium size (bench BIC): ", medbicstats_bench.avgval)

println("Std. Dev rank for medium size (AIC): ", round.(medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(medbicstats.stdval, digits=4))
println("Std. Dev rank for medium size (bench BIC): ", round.(medbicstats_bench.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
println("Freq. Correct for medium size (bench BIC): ", medbicstats_bench.freqcorrect)

println("________________________________________________")

println("BONUS: Statistics if we 'accidentally' choose the incorrect rank")

under_medaicstats = sim_stats(under_medaic21, ranks, sims)
under_medbicstats = sim_stats(under_medbic21, ranks, sims)
under_medbicstats_bench = sim_stats(under_medbic21_bench, ranks, sims)
under_smallaicstats = sim_stats(under_smallaic21, ranks, sims)
under_smallbicstats = sim_stats(under_smallbic21, ranks, sims)
under_smallbicstats_bench = sim_stats(under_smallbic21_bench, ranks, sims)

println("Average rank for small size (AIC): ", under_smallaicstats.avgval)
println("Average rank for small size (BIC): ", under_smallbicstats.avgval)
println("Average rank for small size (bench BIC): ", under_smallbicstats_bench.avgval)

println("Std. Dev rank for small size (AIC): ", round.(under_smallaicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(under_smallbicstats.stdval, digits=4))
println("Std. Dev rank for small size (bench BIC): ", round.(under_smallbicstats_bench.stdval, digits=4))

println("Freq. Correct for small size (AIC): ", under_smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", under_smallbicstats.freqcorrect)
println("Freq. Correct for small size (bench BIC): ", under_smallbicstats_bench.freqcorrect)

println("Average rank for medium size (AIC): ", under_medaicstats.avgval)
println("Average rank for medium size (BIC): ", under_medbicstats.avgval)
println("Average rank for medium size (bench BIC): ", under_medbicstats_bench.avgval)

println("Std. Dev rank for medium size (AIC): ", round.(under_medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(under_medbicstats.stdval, digits=4))
println("Std. Dev rank for medium size (bench BIC): ", round.(under_medbicstats_bench.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", under_medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", under_medbicstats.freqcorrect)
println("Freq. Correct for medium size (bench BIC): ", under_medbicstats_bench.freqcorrect)

println("________________________________________________")

over_medaicstats = sim_stats(over_medaic21, ranks, sims)
over_medbicstats = sim_stats(over_medbic21, ranks, sims)
over_medbicstats_bench = sim_stats(over_medbic21_bench, ranks, sims)
over_smallaicstats = sim_stats(over_smallaic21, ranks, sims)
over_smallbicstats = sim_stats(over_smallbic21, ranks, sims)
over_smallbicstats_bench = sim_stats(over_smallbic21_bench, ranks, sims)

println("Average rank for small size (AIC): ", over_smallaicstats.avgval)
println("Average rank for small size (BIC): ", over_smallbicstats.avgval)
println("Average rank for small size (bench BIC): ", over_smallbicstats_bench.avgval)

println("Std. Dev rank for small size (AIC): ", round.(over_smallaicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(over_smallbicstats.stdval, digits=4))
println("Std. Dev rank for small size (bench BIC): ", round.(over_smallbicstats_bench.stdval, digits=4))

println("Freq. Correct for small size (AIC): ", over_smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", over_smallbicstats.freqcorrect)
println("Freq. Correct for small size (bench BIC): ", over_smallbicstats_bench.freqcorrect)

println("Average rank for medium size (AIC): ", over_medaicstats.avgval)
println("Average rank for medium size (BIC): ", over_medbicstats.avgval)
println("Average rank for medium size (bench BIC): ", over_medbicstats_bench.avgval)

println("Std. Dev rank for medium size (AIC): ", round.(over_medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(over_medbicstats.stdval, digits=4))
println("Std. Dev rank for medium size (bench BIC): ", round.(over_medbicstats_bench.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", over_medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", over_medbicstats.freqcorrect)
println("Freq. Correct for medium size (bench BIC): ", over_medbicstats_bench.freqcorrect)
