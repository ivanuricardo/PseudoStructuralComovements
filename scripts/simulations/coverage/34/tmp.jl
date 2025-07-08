using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250607)

dimvals = [3, 4]
true_rank = [2, 3]
under_rank = [1, 3]
over_rank = [3, 3]

sims = 1000
burnin = 50
obs = 100 + burnin

coef = generate_rrmar_coef(dimvals, true_rank)
delta_true = coef.delta
gamma_true = coef.gamma
u3_true = coef.u3
u4_true = coef.u4

correct_gamma = fill(NaN, 1, sims)
under_gamma = fill(NaN, 1, sims)
over_gamma = fill(NaN, 1, sims)

@showprogress Threads.@threads for i = 1:sims
    data = simulate_rrmar_data(dimvals, true_rank, obs; A=coef, burnin)
    cen_data = data.data .- mean(data.data, dims=2)

    correct_reg = comovement_reg(cen_data, dimvals, true_rank; iters=1000)
    over_reg = comovement_reg(cen_data, dimvals, over_rank; iters=1000)
    under_reg = comovement_reg(cen_data, dimvals, under_rank; iters=1000)

    correct_gamma[:, i] = correct_reg.count
    under_gamma[:, i] = under_reg.count
    over_gamma[:, i] = over_reg.count

end

println("Mean of correct: ", mean(correct_gamma))
println("Mean of over: ", mean(over_gamma))
println("Mean of under: ", mean(under_gamma))

counts = Dict(k => count(==(k), correct_gamma) for k in unique(correct_gamma))

total = length(correct_gamma)
percentages = Dict(k => 100 * v / total for (k, v) in counts)

println("FOR CORRECT GAMMA")
for i in sort(collect(keys(percentages)))
    println("Value $i: $(round(percentages[i]; digits=2))%")
end

################################################################################

counts = Dict(k => count(==(k), under_gamma) for k in unique(under_gamma))

total = length(under_gamma)
percentages = Dict(k => 100 * v / total for (k, v) in counts)

println("FOR UNDER GAMMA")
for i in sort(collect(keys(percentages)))
    println("Value $i: $(round(percentages[i]; digits=2))%")
end

################################################################################

counts = Dict(k => count(==(k), over_gamma) for k in unique(over_gamma))

total = length(over_gamma)
percentages = Dict(k => 100 * v / total for (k, v) in counts)

println("FOR OVER GAMMA")
for i in sort(collect(keys(percentages)))
    println("Value $i: $(round(percentages[i]; digits=2))%")
end
