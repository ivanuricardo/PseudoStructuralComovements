using DrWatson
@quickactivate :PseudoStructuralComovements

ranks = [3, 1]
sims = 100
sim_results1 = load(datadir("threebyfour/31_results.jld2"))
sim_results2 = load(datadir("threebyfour/under_31_results.jld2"))
sim_results3 = load(datadir("threebyfour/second_31_results.jld2"))
smallaic = sim_results1["smallaic"]
smallbic = sim_results1["smallbic"]
medaic = sim_results1["medaic"]
medbic = sim_results1["medbic"]

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

################################################################################

println("________________________________________________")

println("BONUS: Statistics if we 'accidentally' choose the incorrect rank")

under_smallaic = sim_results2["under_smallaic"]
under_smallbic = sim_results2["under_smallbic"]
medaic = sim_results1["medaic"]
medbic = sim_results1["medbic"]

medaicstats = sim_stats(medaic, ranks, sims)
medbicstats = sim_stats(medbic, ranks, sims)
under_smallaicstats = sim_stats(under_smallaic, ranks, sims)
under_smallbicstats = sim_stats(under_smallbic, ranks, sims)

println("Average rank for under_small size (AIC): ", under_smallaicstats.avgval)
println("Average rank for under_small size (BIC): ", under_smallbicstats.avgval)

println("Std. Dev rank for under_small size (AIC): ", round.(under_smallaicstats.stdval, digits=4))
println("Std. Dev rank for under_small size (BIC): ", round.(under_smallbicstats.stdval, digits=4))

println("Freq. Correct for under_small size (AIC): ", under_smallaicstats.freqcorrect)
println("Freq. Correct for under_small size (BIC): ", under_smallbicstats.freqcorrect)

println("Average rank for medium size (AIC): ", medaicstats.avgval)
println("Average rank for medium size (BIC): ", medbicstats.avgval)

println("Std. Dev rank for medium size (AIC): ", round.(medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(medbicstats.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)

println("________________________________________________")

second_smallaic = sim_results3["second_smallaic"]
second_smallbic = sim_results3["second_smallbic"]
medaic = sim_results1["medaic"]
medbic = sim_results1["medbic"]

medaicstats = sim_stats(medaic, ranks, sims)
medbicstats = sim_stats(medbic, ranks, sims)
second_smallaicstats = sim_stats(second_smallaic, ranks, sims)
second_smallbicstats = sim_stats(second_smallbic, ranks, sims)

println("Average rank for second_small size (AIC): ", second_smallaicstats.avgval)
println("Average rank for second_small size (BIC): ", second_smallbicstats.avgval)

println("Std. Dev rank for second_small size (AIC): ", round.(second_smallaicstats.stdval, digits=4))
println("Std. Dev rank for second_small size (BIC): ", round.(second_smallbicstats.stdval, digits=4))

println("Freq. Correct for second_small size (AIC): ", second_smallaicstats.freqcorrect)
println("Freq. Correct for second_small size (BIC): ", second_smallbicstats.freqcorrect)

println("Average rank for medium size (AIC): ", medaicstats.avgval)
println("Average rank for medium size (BIC): ", medbicstats.avgval)

println("Std. Dev rank for medium size (AIC): ", round.(medaicstats.stdval, digits=4))
println("Std. Dev rank for medium size (BIC): ", round.(medbicstats.stdval, digits=4))

println("Freq. Correct for medium size (AIC): ", medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
