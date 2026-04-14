using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250723)

dimvals = [3, 4]
true_rank = [2, 3]
under_rank = [1, 3]
over_rank = [3, 3]

sims = 1000
burnin = 50
obs = 250 + burnin

coef = generate_rrmar_coef(dimvals, true_rank)
delta_true = coef.delta
gamma_true = coef.gamma
u3_true = coef.u3
u4_true = coef.u4

correct_gamma = fill(NaN, 3, sims)
under_gamma = fill(NaN, 3, sims)
over_gamma = fill(NaN, 3, sims)
ps_ll = fill(NaN, 1, sims)
comove_iters = fill(NaN, 1, sims)

correct_rrmar = fill(NaN, 3, sims)
under_rrmar = fill(NaN, 3, sims)
over_rrmar = fill(NaN, 3, sims)
rrmar_ll = fill(NaN, 1, sims)

@showprogress Threads.@threads for i = 1:sims
    data = simulate_rrmar_data(dimvals, true_rank, obs; A=coef, burnin, matrix_err=true)

    correct_reg = comovement_reg(data.data, dimvals, true_rank; iters=1000, num_starts=100, num_selected=10)
    over_reg = comovement_reg(data.data, dimvals, over_rank; iters=1000, num_starts=100, num_selected=10)
    under_reg = comovement_reg(data.data, dimvals, under_rank; iters=1000, num_starts=100, num_selected=10)

    rrmar_correct = rrmar(data.data, dimvals, true_rank)
    rrmar_over = rrmar(data.data, dimvals, over_rank)
    rrmar_under = rrmar(data.data, dimvals, under_rank)

    correct_gamma[:, i] = correct_reg.gamma_est[2:end]
    under_gamma[:, i] = under_reg.gamma_est[2:end]
    over_gamma[:, i] = over_reg.gamma_est[2:end]
    ps_ll[i] = correct_reg.res.minimum
    comove_iters[i] = correct_reg.res.iterations

    correct_rrmar[:, i] = rrmar_correct.gamma_est[2:end]
    under_rrmar[:, i] = rrmar_under.gamma_est[2:end]
    over_rrmar[:, i] = rrmar_over.gamma_est[2:end]
    rrmar_ll[i] = rrmar_correct.ll

end

save(datadir("coverage/34/gamma_comparison_results250.jld2"), Dict(
    "correct_gamma" => correct_gamma,
    "under_gamma" => under_gamma,
    "over_gamma" => over_gamma,
    "correct_rrmar" => correct_rrmar,
    "under_rrmar" => under_rrmar,
    "over_rrmar" => over_rrmar,
    "gamma_true" => gamma_true,
    "ps_ll" => ps_ll,
    "rrmar_ll" => rrmar_ll,
    "comove_iters" => comove_iters,
))

h1 = StatsPlots.density(
    correct_rrmar[1, :];
    legend=false,
    ylabel="Density",
    linewidth=3,
    yguidefont=16   # <-- increase y-axis label font size
)
StatsPlots.density!(under_rrmar[1, :]; linewidth=3)
StatsPlots.density!(over_rrmar[1, :]; linewidth=3)
vline!([gamma_true[2]]; linewidth=3)

h2 = StatsPlots.density(
    correct_gamma[1, :];
    legend=false,
    ylabel="Density",
    linewidth=3,
    yguidefont=16   # <-- increase y-axis label font size
)
StatsPlots.density!(under_gamma[1, :]; linewidth=3)
StatsPlots.density!(over_gamma[1, :]; linewidth=3)
vline!([gamma_true[2]]; linewidth=3)
