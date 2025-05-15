using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250421)

dimvals = [3, 4]
true_ranks = [2, 2]
under_rank = [2, 1]
over_rank = [2, 3]

sims = 100
burnin = 50
obs = 100 + burnin

coef = generate_rrmar_coef(dimvals, true_ranks)
u1_true = coef.u1
u2_true = coef.u2
u3_true = coef.u3
u4_true = coef.u4

delta_true = nullspace(u1_true') ./ nullspace(u1_true')[1]
gamma_true = nullspace(u2_true') * inv(nullspace(u2_true')[1:2, 1:2])
u3_rot = u3_true * inv(u3_true[1:2, 1:2])
u4_rot = u4_true * inv(u4_true[1:2, 1:2])

correct_delta = fill(NaN, 2, sims)
under_delta = fill(NaN, 2, sims)
over_delta = fill(NaN, 2, sims)

correct_cov = fill(NaN, 2, sims)
under_cov = fill(NaN, 2, sims)
over_cov = fill(NaN, 2, sims)

@showprogress for i = 1:sims
    data = simulate_rrmar_data(dimvals, true_ranks, obs; A=coef, burnin)
    cen_data = data.data .- mean(data.data, dims=2)

    @time correct_reg = comovement_reg(cen_data, dimvals, true_ranks, iters=100)
    println(correct_reg.res.iterations)
    if correct_reg.res.g_residual > 1.0
        println("Not converged")
    end
    over_reg = comovement_reg(cen_data, dimvals, over_rank, iters=100)
    under_reg = comovement_reg(cen_data, dimvals, under_rank, iters=100)

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

save(datadir("coverage/delta_cov_results100.jld2"), Dict(
    "correct_delta" => correct_delta,
    "under_delta" => under_delta,
    "over_delta" => over_delta,
    "correct_cov" => correct_cov,
    "under_cov" => under_cov,
    "over_cov" => over_cov,
    "delta_true" => delta_true
))

