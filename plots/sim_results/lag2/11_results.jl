using DrWatson
@quickactivate :PseudoStructuralComovements

ranks = [1, 1, 2]
sims = 100
sim_results = load(datadir("lag2/11_results.jld2"))
smallaic = sim_results["smallaic"]
smallbic = sim_results["smallbic"]
medaic = sim_results["medaic"]
medbic = sim_results["medbic"]

medaicstats = sim_stats(medaic, ranks, sims)
medbicstats = sim_stats(medbic, ranks, sims)
smallaicstats = sim_stats(smallaic, ranks, sims)
smallbicstats = sim_stats(smallbic, ranks, sims)

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
