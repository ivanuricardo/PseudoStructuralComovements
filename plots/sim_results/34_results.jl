using DrWatson
@quickactivate :PseudoStructuralComovements

ranks = [3, 4]
sims = 100
sim_results = load(datadir("threebyfour/34_results.jld2"))
smallaic = sim_results["smallaic"]
smallbic = sim_results["smallbic"]
smallbic_bench = sim_results["smallbic_bench"]
medbic_bench = sim_results["medbic_bench"]
medaic = sim_results["medaic"]
medbic = sim_results["medbic"]

medaicstats = sim_stats(medaic, ranks, sims)
medbicstats = sim_stats(medbic, ranks, sims)
medbicstats_bench = sim_stats(medbic_bench, ranks, sims)
smallaicstats = sim_stats(smallaic, ranks, sims)
smallbicstats = sim_stats(smallbic, ranks, sims)
smallbicstats_bench = sim_stats(smallbic_bench, ranks, sims)

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

