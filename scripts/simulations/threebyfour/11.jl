using DrWatson
@quickactivate :PseudoStructuralComovements
using RCall
R"""
source("r_helpers.R")
"""
Random.seed!(20250713)
const R_LOCK = ReentrantLock()

sims = 100
dimvals = [3, 4]
ranks = [1, 1]
snr = 0.7

smallaic11 = fill(NaN, 2, sims)
smallbic11 = fill(NaN, 2, sims)
medaic11 = fill(NaN, 2, sims)
medbic11 = fill(NaN, 2, sims)

smallbic11_bench = fill(NaN, 2, sims)
medbic11_bench = fill(NaN, 2, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin, matrix_err=true)
    med_bench_data = reshape(medmar.data', (medobs, dimvals[1], dimvals[2]))
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin, matrix_err=true)
    small_bench_data = reshape(smallmar.data', (smallobs, dimvals[1], dimvals[2]))

    smallicest = rank_selection(smallmar.data, dimvals; iters=1000)
    smallaic11[:, s] .= smallicest.aic_sel[1:2]
    smallbic11[:, s] .= smallicest.bic_sel[1:2]

    medicest = rank_selection(medmar.data, dimvals; iters=1000)
    medaic11[:, s] .= medicest.aic_sel[1:2]
    medbic11[:, s] .= medicest.bic_sel[1:2]

    lock(R_LOCK) do
        small_bench = R"""
        d1 = $dimvals[1]
        d2 = $dimvals[2]
        small_data <- $small_bench_data
        small_selected_rank <- r_rank_selection(small_data, d1, d2)
        """
        @rget small_selected_rank
        smallbic11_bench[:, s] .= small_selected_rank[:selected_ranks]
    end

    lock(R_LOCK) do
        med_bench = R"""
        d1 = $dimvals[1]
        d2 = $dimvals[2]
        med_data <- $med_bench_data
        med_selected_rank <- r_rank_selection(med_data, d1, d2)
        """
        @rget med_selected_rank
        medbic11_bench[:, s] .= med_selected_rank[:selected_ranks]
    end
end

save(datadir("threebyfour/11_results.jld2"), Dict(
    "smallaic" => smallaic11,
    "smallbic" => smallbic11,
    "medaic" => medaic11,
    "medbic" => medbic11,
    "smallbic_bench" => smallbic11_bench,
    "medbic_bench" => medbic11_bench,
))

medaicstats = sim_stats(medaic11, ranks, sims)
medbicstats = sim_stats(medbic11, ranks, sims)
smallaicstats = sim_stats(smallaic11, ranks, sims)
smallbicstats = sim_stats(smallbic11, ranks, sims)
smallbicstats_bench = sim_stats(smallbic11_bench, ranks, sims)
medbicstats_bench = sim_stats(medbic11_bench, ranks, sims)

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
