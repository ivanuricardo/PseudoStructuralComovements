using DrWatson
@quickactivate :PseudoStructuralComovements
using RCall
Random.seed!(20250723)

sims = 100
dimvals = [3, 4]
ranks = [1, 4]
snr = 0.7
p = 2

smallaic14 = fill(NaN, 3, sims)
smallbic14 = fill(NaN, 3, sims)
medaic14 = fill(NaN, 3, sims)
medbic14 = fill(NaN, 3, sims)

burnin = 100
smallobs = 100
medobs = 250

A = generate_rrmar_coef(dimvals, ranks; p=2)

@showprogress Threads.@threads for s = 1:sims
    medmar = simulate_rrmar_data(dimvals, ranks, medobs + burnin; A, snr, burnin, matrix_err=true, p=2)
    smallmar = simulate_rrmar_data(dimvals, ranks, smallobs + burnin; A, snr, burnin, matrix_err=true, p=2)

    smallicest = rank_selection(smallmar.data, dimvals; iters=1000, pmax=2)
    smallaic14[:, s] .= smallicest.aic_sel[1:3]
    smallbic14[:, s] .= smallicest.bic_sel[1:3]

    medicest = rank_selection(medmar.data, dimvals; iters=1000, pmax=2)
    medaic14[:, s] .= medicest.aic_sel[1:3]
    medbic14[:, s] .= medicest.bic_sel[1:3]
end

save(datadir("lag1/14_results.jld2"), Dict(
    "smallaic" => smallaic14,
    "smallbic" => smallbic14,
    "medaic" => medaic14,
    "medbic" => medbic14,
))

push!(ranks, 1)
medaicstats = sim_stats(medaic14, ranks, sims)
medbicstats = sim_stats(medbic14, ranks, sims)
smallaicstats = sim_stats(smallaic14, ranks, sims)
smallbicstats = sim_stats(smallbic14, ranks, sims)

println("Average rank for small size (AIC): ", smallaicstats.avgval)
println("Average lag for small size (AIC): ", smallaicstats.avgval)
println("Average rank for small size (BIC): ", smallbicstats.avgval)
println("Average lag for small size (BIC): ", smallbicstats.avgval)

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
