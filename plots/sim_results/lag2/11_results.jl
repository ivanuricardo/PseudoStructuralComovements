using DrWatson
@quickactivate :PseudoStructuralComovements

ranks = [1, 1, 2]
sims = 100
sim_results = load(datadir("lag2/11_results.jld2"))
smallaic = sim_results["smallaic"]
smallbic = sim_results["smallbic"]
smallebic = sim_results["smallebic"]
medaic = sim_results["medaic"]
medbic = sim_results["medbic"]
medebic = sim_results["medebic"]

medaicstats = sim_stats(medaic, ranks, sims)
medbicstats = sim_stats(medbic, ranks, sims)
medebicstats = sim_stats(medebic, ranks, sims)
smallaicstats = sim_stats(smallaic, ranks, sims)
smallbicstats = sim_stats(smallbic, ranks, sims)
smallebicstats = sim_stats(smallebic, ranks, sims)

println("MAD for small size (AIC): ", smallaicstats.mad)
println("MAD for small size (BIC): ", smallbicstats.mad)
println("MAD for small size (EBIC): ", smallebicstats.mad)

println("Freq. Correct for small size (AIC): ", smallaicstats.freqcorrect)
println("Freq. Correct for small size (BIC): ", smallbicstats.freqcorrect)
println("Freq. Correct for small size (EBIC): ", smallebicstats.freqcorrect)

println("MAD for medium size (AIC): ", medaicstats.mad)
println("MAD for medium size (BIC): ", medbicstats.mad)
println("MAD for medium size (EBIC): ", medebicstats.mad)

println("Freq. Correct for medium size (AIC): ", medaicstats.freqcorrect)
println("Freq. Correct for medium size (BIC): ", medbicstats.freqcorrect)
println("Freq. Correct for medium size (EBIC): ", medebicstats.freqcorrect)
