using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250530)

dimvals = [3, 4]
true_ranks = [2, 2]
under_rank = [2, 1]
over_rank = [2, 3]

sims = 250
burnin = 50
obs = 100 + burnin

coef = generate_rrmar_coef(dimvals, true_ranks)
delta_true = coef.delta
gamma_true = coef.gamma
u3_true = coef.u3
u4_true = coef.u4

correct_delta = fill(NaN, 2, sims)
under_delta = fill(NaN, 2, sims)
over_delta = fill(NaN, 2, sims)

correct_cov = fill(NaN, 2, sims)
under_cov = fill(NaN, 2, sims)
over_cov = fill(NaN, 2, sims)

@showprogress Threads.@threads for i = 1:sims
    data = simulate_rrmar_data(dimvals, true_ranks, obs; A=coef, burnin)
    cen_data = data.data .- mean(data.data, dims=2)

    correct_reg = comovement_reg(cen_data, dimvals, true_ranks; iters=100)
    over_reg = comovement_reg(cen_data, dimvals, over_rank; iters=100)
    under_reg = comovement_reg(cen_data, dimvals, under_rank; iters=100)

    correct_delta[:, i] = correct_reg.delta_est[2:end]
    under_delta[:, i] = under_reg.delta_est[2:end]
    over_delta[:, i] = over_reg.delta_est[2:end]

    correct_upper = correct_reg.delta_est[2:end] + 1.96 .* correct_reg.delta_stderr
    correct_lower = correct_reg.delta_est[2:end] - 1.96 .* correct_reg.delta_stderr
    correct_cov[:, i] = correct_lower .< delta_true[2:end] .< correct_upper

    under_upper = under_reg.delta_est[2:end] + 1.96 .* under_reg.delta_stderr
    under_lower = under_reg.delta_est[2:end] - 1.96 .* under_reg.delta_stderr
    under_cov[:, i] = under_lower .< delta_true[2:end] .< under_upper

    over_upper = over_reg.delta_est[2:end] + 1.96 .* over_reg.delta_stderr
    over_lower = over_reg.delta_est[2:end] - 1.96 .* over_reg.delta_stderr
    over_cov[:, i] = over_lower .< delta_true[2:end] .< over_upper
end

density(filter(!isnan, correct_delta[2, :]))
density(filter(!isnan, over_delta[1, :]))
density(filter(!isnan, under_delta[1, :]))

save(datadir("coverage/delta_cov_results100.jld2"), Dict(
    "correct_delta" => correct_delta,
    "under_delta" => under_delta,
    "over_delta" => over_delta,
    "correct_cov" => correct_cov,
    "under_cov" => under_cov,
    "over_cov" => over_cov,
    "delta_true" => delta_true
))

