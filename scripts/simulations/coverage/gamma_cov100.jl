using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250421)

dimvals = [3, 4]
true_ranks = [2, 3]
under_rank = [1, 3]
over_rank = [3, 3]

sims = 1000
burnin = 50
obs = 100 + burnin

coef = generate_rrmar_coef(dimvals, true_ranks)
u1_true = coef.u1
u2_true = coef.u2
u3_true = coef.u3
u4_true = coef.u4

delta_true = nullspace(u1_true') ./ nullspace(u1_true')[1]
gamma_true = nullspace(u2_true') ./ nullspace(u2_true')[1]
u3_rot = u3_true * inv(u3_true[1:2, 1:2])
u4_rot = u4_true * inv(u4_true[1:3, 1:3])

correct_gamma = fill(NaN, 3, sims)
under_gamma = fill(NaN, 3, sims)
over_gamma = fill(NaN, 3, sims)

correct_cov = fill(NaN, 3, sims)
under_cov = fill(NaN, 3, sims)
over_cov = fill(NaN, 3, sims)

#=@showprogress Threads.@threads for i = 1:sims=#
@showprogress for i = 1:sims
    data = simulate_rrmar_data(dimvals, true_ranks, obs; A=coef, burnin)
    cen_data = data.data .- mean(data.data, dims=2)

    correct_reg = comovement_reg(cen_data, dimvals, true_ranks, iters=500)
    over_reg = comovement_reg(cen_data, dimvals, over_rank, iters=500)
    under_reg = comovement_reg(cen_data, dimvals, under_rank, iters=500)

    correct_gamma[:, i] = correct_reg.gamma_est[2:end]
    under_gamma[:, i] = under_reg.gamma_est[2:end]
    over_gamma[:, i] = over_reg.gamma_est[2:end]

    correct_upper = correct_reg.gamma_est[2:end] + 1.96 .* correct_reg.gamma_stderr
    correct_lower = correct_reg.gamma_est[2:end] - 1.96 .* correct_reg.gamma_stderr
    correct_cov[:, i] = correct_lower .< gamma_true[2:end] .< correct_upper

    under_upper = under_reg.gamma_est[2:end] + 1.96 .* under_reg.gamma_stderr
    under_lower = under_reg.gamma_est[2:end] - 1.96 .* under_reg.gamma_stderr
    under_cov[:, i] = under_lower .< gamma_true[2:end] .< under_upper

    over_upper = over_reg.gamma_est[2:end] + 1.96 .* over_reg.gamma_stderr
    over_lower = over_reg.gamma_est[2:end] - 1.96 .* over_reg.gamma_stderr
    over_cov[:, i] = over_lower .< gamma_true[2:end] .< over_upper
end

save(datadir("coverage/gamma_cov_results100.jld2"), Dict(
    "correct_gamma" => correct_gamma,
    "under_gamma" => under_gamma,
    "over_gamma" => over_gamma,
    "correct_cov" => correct_cov,
    "under_cov" => under_cov,
    "over_cov" => over_cov,
    "gamma_true" => gamma_true
))

h1 = density(true_gamma[1, :], label="Correct rank")
density!(under_gamma[1, :], label="Underestimated rank")
density!(over_gamma[1, :], label="Overestimated rank")
vline!([gamma_true[2]], label="True value")
vline!([mean(true_gamma[1, :])], label="Correct rank mean")
vline!([mean(under_gamma[1, :])], label="Under rank mean")
vline!([mean(over_gamma[1, :])], label="Over rank mean")

h2 = density(true_gamma[2, :], label="Correct rank")
density!(under_gamma[2, :], label="Underestimated rank")
density!(over_gamma[2, :], label="Overestimated rank")
vline!([gamma_true[3]], label="True value")

h3 = density(true_gamma[3, :], label="Correct rank")
density!(under_gamma[3, :], label="Underestimated rank")
density!(over_gamma[3, :], label="Overestimated rank")
vline!([gamma_true[4]], label="True value")
