using DrWatson
@quickactivate :PseudoStructuralComovements

ranks = [1, 4]
sims = 100
sim_results = load(datadir("threebyfour/14_bench_results.jld2"))
smallbic = sim_results["smallbic"]
medbic = sim_results["medbic"]

medbicstats = sim_stats(medbic, ranks, sims)
smallbicstats = sim_stats(smallbic, ranks, sims)

println("Average rank for small size (BIC): ", smallbicstats.avgval)
println("Std. Dev rank for small size (BIC): ", round.(smallbicstats.stdval, digits=4))
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)
println("Average rank for medium size (BIC): ", medbicstats.avgval)
println("Std. Dev rank for medium size (BIC): ", round.(medbicstats.stdval, digits=4))
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
