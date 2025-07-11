using DrWatson
@quickactivate :PseudoStructuralComovements

ranks = [2, 1]
sims = 100
sim_results = load(datadir("threebyfour/21_results.jld2"))
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

println("________________________________________________")

println("BONUS: Statistics if we 'accidentally' choose the incorrect rank")

under_sim_results = load(datadir("threebyfour/under_21_results.jld2"))
under_smallaic = under_sim_results["under_smallaic"]
under_smallbic = under_sim_results["under_smallbic"]
under_smallbic_bench = under_sim_results["under_smallbic_bench"]
under_medbic_bench = under_sim_results["under_medbic_bench"]
under_medaic = under_sim_results["under_medaic"]
under_medbic = under_sim_results["under_medbic"]

under_medaicstats = sim_stats(under_medaic, ranks, sims)
under_medbicstats = sim_stats(under_medbic, ranks, sims)
under_medbicstats_bench = sim_stats(under_medbic_bench, ranks, sims)
under_smallaicstats = sim_stats(under_smallaic, ranks, sims)
under_smallbicstats = sim_stats(under_smallbic, ranks, sims)
under_smallbicstats_bench = sim_stats(under_smallbic_bench, ranks, sims)

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

over_sim_results = load(datadir("threebyfour/over_21_results.jld2"))
over_smallaic = over_sim_results["over_smallaic"]
over_smallbic = over_sim_results["over_smallbic"]
over_smallbic_bench = over_sim_results["over_smallbic_bench"]
over_medbic_bench = over_sim_results["over_medbic_bench"]
over_medaic = over_sim_results["over_medaic"]
over_medbic = over_sim_results["over_medbic"]

over_medaicstats = sim_stats(over_medaic, ranks, sims)
over_medbicstats = sim_stats(over_medbic, ranks, sims)
over_medbicstats_bench = sim_stats(over_medbic_bench, ranks, sims)
over_smallaicstats = sim_stats(over_smallaic, ranks, sims)
over_smallbicstats = sim_stats(over_smallbic, ranks, sims)
over_smallbicstats_bench = sim_stats(over_smallbic_bench, ranks, sims)

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
