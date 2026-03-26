using DrWatson
@quickactivate :PseudoStructuralComovements

ranks = [2, 1]
sims = 1000
sim_results = load(datadir("threebyfour/21_results.jld2"))
smallaic = sim_results["smallaic"]
smallbic = sim_results["smallbic"]
smallebic = sim_results["smallebic"]
smallbench = sim_results["smallbench"]
medaic = sim_results["medaic"]
medbic = sim_results["medbic"]
medebic = sim_results["medebic"]
medbench = sim_results["medbench"]

medaicstats = sim_stats(medaic, ranks, sims)
medbicstats = sim_stats(medbic, ranks, sims)
medebicstats = sim_stats(medebic, ranks, sims)
medbenchstats = sim_stats(medbench, ranks, sims)

smallaicstats = sim_stats(smallaic, ranks, sims)
smallbicstats = sim_stats(smallbic, ranks, sims)
smallebicstats = sim_stats(smallebic, ranks, sims)
smallbenchstats = sim_stats(smallbench, ranks, sims)

println("MAD for small size (AIC): ", smallaicstats.mad)
println("MAD for small size (BIC): ", smallbicstats.mad)
println("MAD for small size (EBIC): ", smallebicstats.mad)
println("MAD for small size (Bench EBIC): ", smallbenchstats.mad)

println("Freq. Correct for small size (AIC): ", smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)
println("Freq. Correct for small size (EBIC): ", smallebicstats.freqcorrect)
println("Freq. Correct for small size (Bench EBIC): ", smallbenchstats.freqcorrect)

println("MAD for medium size (AIC): ", medaicstats.mad)
println("MAD for medium size (BIC): ", medbicstats.mad)
println("MAD for medium size (EBIC): ", medebicstats.mad)
println("MAD for medium size (Bench EBIC): ", medbenchstats.mad)

println("Freq. Correct for medium size (AIC): ", medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
println("Freq. Correct for medium size (EBIC): ", medebicstats.freqcorrect)
println("Freq. Correct for medium size (Bench EBIC): ", medbenchstats.freqcorrect)

println("________________________________________________")

println("BONUS: Statistics if we 'accidentally' choose the incorrect rank")

under_sim_results = load(datadir("threebyfour/under_21_results.jld2"))
under_smallaic = under_sim_results["under_smallaic"]
under_smallbic = under_sim_results["under_smallbic"]
under_smallebic = under_sim_results["under_smallebic"]
under_smallbench = under_sim_results["under_smallbench"]
under_medaic = under_sim_results["under_medaic"]
under_medbic = under_sim_results["under_medbic"]
under_medebic = under_sim_results["under_medebic"]
under_medbench = under_sim_results["under_medbench"]

under_medaicstats = sim_stats(under_medaic, ranks, sims)
under_medbicstats = sim_stats(under_medbic, ranks, sims)
under_medebicstats = sim_stats(under_medebic, ranks, sims)
under_medbenchstats = sim_stats(under_medbench, ranks, sims)

under_smallaicstats = sim_stats(under_smallaic, ranks, sims)
under_smallbicstats = sim_stats(under_smallbic, ranks, sims)
under_smallebicstats = sim_stats(under_smallebic, ranks, sims)
under_smallbenchstats = sim_stats(under_smallbench, ranks, sims)

println("MAD for small size (AIC): ", under_smallaicstats.mad)
println("MAD for small size (BIC): ", under_smallbicstats.mad)
println("MAD for small size (EBIC): ", under_smallebicstats.mad)
println("MAD for small size (Bench EBIC): ", under_smallbenchstats.mad)

println("Freq. Correct for small size (AIC): ", under_smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", under_smallbicstats.freqcorrect)
println("Freq. Correct for small size (EBIC): ", under_smallebicstats.freqcorrect)
println("Freq. Correct for small size (Bench EBIC): ", under_smallbenchstats.freqcorrect)

println("MAD for medium size (AIC): ", under_medaicstats.mad)
println("MAD for medium size (BIC): ", under_medbicstats.mad)
println("MAD for medium size (EBIC): ", under_medebicstats.mad)
println("MAD for medium size (Bench EBIC): ", under_medbenchstats.mad)

println("Freq. Correct for medium size (AIC): ", under_medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", under_medbicstats.freqcorrect)
println("Freq. Correct for medium size (EBIC): ", under_medebicstats.freqcorrect)
println("Freq. Correct for medium size (Bench EBIC): ", under_medbenchstats.freqcorrect)

println("________________________________________________")

over_sim_results = load(datadir("threebyfour/over_21_results.jld2"))
over_smallaic = over_sim_results["over_smallaic"]
over_smallbic = over_sim_results["over_smallbic"]
over_smallebic = over_sim_results["over_smallebic"]
over_smallbench = over_sim_results["over_smallbench"]
over_medaic = over_sim_results["over_medaic"]
over_medbic = over_sim_results["over_medbic"]
over_medebic = over_sim_results["over_medebic"]
over_medbench = over_sim_results["over_medbench"]

over_medaicstats = sim_stats(over_medaic, ranks, sims)
over_medbicstats = sim_stats(over_medbic, ranks, sims)
over_medebicstats = sim_stats(over_medebic, ranks, sims)
over_medbenchstats = sim_stats(over_medbench, ranks, sims)

over_smallaicstats = sim_stats(over_smallaic, ranks, sims)
over_smallbicstats = sim_stats(over_smallbic, ranks, sims)
over_smallebicstats = sim_stats(over_smallebic, ranks, sims)
over_smallbenchstats = sim_stats(over_smallbench, ranks, sims)

println("MAD for small size (AIC): ", over_smallaicstats.mad)
println("MAD for small size (BIC): ", over_smallbicstats.mad)
println("MAD for small size (EBIC): ", over_smallebicstats.mad)
println("MAD for small size (Bench EBIC): ", over_smallbenchstats.mad)

println("Freq. Correct for small size (AIC): ", over_smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", over_smallbicstats.freqcorrect)
println("Freq. Correct for small size (EBIC): ", over_smallebicstats.freqcorrect)
println("Freq. Correct for small size (Bench EBIC): ", over_smallbenchstats.freqcorrect)

println("MAD for medium size (AIC): ", over_medaicstats.mad)
println("MAD for medium size (BIC): ", over_medbicstats.mad)
println("MAD for medium size (EBIC): ", over_medebicstats.mad)
println("MAD for medium size (Bench EBIC): ", over_medbenchstats.mad)

println("Freq. Correct for medium size (AIC): ", over_medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", over_medbicstats.freqcorrect)
println("Freq. Correct for medium size (EBIC): ", over_medebicstats.freqcorrect)
println("Freq. Correct for medium size (Bench EBIC): ", over_medbenchstats.freqcorrect)
