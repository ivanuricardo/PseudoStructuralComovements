using DrWatson
@quickactivate :PseudoStructuralComovements
using RCall
Random.seed!(20250729)

function sim_results(aic, bic, ebic, ranks, sims)

    aic = aic[:, .!vec(all(isnan.(aic), dims=1))]
    bic = bic[:, .!vec(all(isnan.(bic), dims=1))]
    ebic = ebic[:, .!vec(all(isnan.(ebic), dims=1))]

    aicstats = sim_stats(aic, ranks, sims)
    bicstats = sim_stats(bic, ranks, sims)
    ebicstats = sim_stats(ebic, ranks, sims)

    println("MAD for (AIC): ", aicstats.mad)
    println("MAD for (BIC): ", bicstats.mad)
    println("MAD for (EBIC): ", ebicstats.mad)

    println("Std. Dev rank for (AIC): ", round.(aicstats.stdval, digits=4))
    println("Std. Dev rank for (BIC): ", round.(bicstats.stdval, digits=4))
    println("Std. Dev rank for (EBIC): ", round.(bicstats.stdval, digits=4))

    println("Freq. Correct for (AIC): ", aicstats.freqcorrect)
    println("Freq. Correct for (BIC): ", bicstats.freqcorrect)
    println("Freq. Correct for (EBIC): ", bicstats.freqcorrect)

end

sims = 100
dimvals = [9, 4]
ranks = [2, 1]
snr = 0.7

aic_val = fill(NaN, 2, sims)
bic_val = fill(NaN, 2, sims)
ebic_val = fill(NaN, 2, sims)

burnin = 100
obs = 180

A = generate_rrmar_coef(dimvals, ranks)
A.sorted_eigs

@showprogress Threads.@threads for s = 1:sims
    mar = simulate_rrmar_data(dimvals, ranks, obs + burnin; A, snr, burnin, matrix_err=true)

    icest = rank_selection(mar.data, dimvals; num_starts=20, num_selected=3)
    # icest2 = rrmar_ic(mar.data, dimvals)
    aic_val[:, s] .= icest.aic_sel[1:2]
    bic_val[:, s] .= icest.bic_sel[1:2]
    ebic_val[:, s] .= icest.ebic_sel[1:2]
    if s % 25 == 0
        println("After $s simulation runs:")
        println("______________________________________")
        sim_results(aic_val, bic_val, ebic_val, ranks, s)
        println("______________________________________")
    end

end

sim_results(aic_val, bic_val, ebic_val, ranks, sims)
