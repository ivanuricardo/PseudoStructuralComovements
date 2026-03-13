using DrWatson
@quickactivate :PseudoStructuralComovements
using RCall
Random.seed!(20250723)

sims = 100
dimvals = [9, 4]
ranks = [2, 1]
snr = 0.7

smallaic11 = fill(NaN, 3, sims)
smallbic11 = fill(NaN, 3, sims)
smallhqc11 = fill(NaN, 3, sims)

burnin = 100
smallobs = 180

A = generate_rrmar_coef(dimvals, ranks)

@showprogress Threads.@threads for s = 1:sims
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin, matrix_err=true)

    smallicest = rank_selection(smallmar.data, dimvals; iters=2, pmax=1)
    smallaic11[:, s] .= smallicest.aic_sel[1:3]
    smallbic11[:, s] .= smallicest.bic_sel[1:3]
    smallhqc11[:, s] .= smallicest.hqc_sel[1:3]
end

save(datadir("lag1/94_results.jld2"), Dict(
    "smallaic" => smallaic11,
    "smallbic" => smallbic11,
    "smallhqc" => smallhqc11,
))

push!(ranks, 1)
smallaicstats = sim_stats(smallaic11, ranks, sims)
smallbicstats = sim_stats(smallbic11, ranks, sims)
smallhqcstats = sim_stats(smallhqc11, ranks, sims)

println("MAD for small size (AIC): ", smallaicstats.mad)
println("MAD for small size (BIC): ", smallbicstats.mad)
println("MAD for small size (HQC): ", smallhqcstats.mad)

println("Std. Dev rank for small size (AIC): ", round.(smallaicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(smallbicstats.stdval, digits=4))
println("Std. Dev rank for small size (BIC): ", round.(smallhqcstats.stdval, digits=4))

println("Freq. Correct for small size (AIC): ", smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)
println("Freq. Correct for small size (HQC): ", smallhqcstats.freqcorrect)

