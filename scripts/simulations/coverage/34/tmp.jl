using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250607)

dimvals = [3, 4]
true_rank = [2, 4]
under_rank = [2, 1]

sims = 500
burnin = 50
obs = 250 + burnin

coef = generate_rrmar_coef(dimvals, true_rank)
delta_true = coef.delta
gamma_true = coef.gamma
u3_true = coef.u3
u4_true = coef.u4

correct_delta = fill(NaN, 2, sims)
under_delta = fill(NaN, 2, sims)

correct_cov = fill(NaN, 2, sims)
under_cov = fill(NaN, 2, sims)

@showprogress Threads.@threads for i = 1:sims
    data = simulate_rrmar_data(dimvals, true_rank, obs; A=coef, burnin, snr=0.3)
    cen_data = data.data .- mean(data.data, dims=2)

    correct_reg = comovement_reg(cen_data, dimvals, true_rank; iters=1000)
    under_reg = comovement_reg(cen_data, dimvals, under_rank; iters=1000)

    correct_delta[:, i] = correct_reg.delta_est[2:end]
    under_delta[:, i] = under_reg.delta_est[2:end]

    correct_upper = correct_reg.delta_est[2:end] + 1.96 .* correct_reg.delta_stderr
    correct_lower = correct_reg.delta_est[2:end] - 1.96 .* correct_reg.delta_stderr
    correct_cov[:, i] = correct_lower .< delta_true[2:end] .< correct_upper

    under_upper = under_reg.delta_est[2:end] + 1.96 .* under_reg.delta_stderr
    under_lower = under_reg.delta_est[2:end] - 1.96 .* under_reg.delta_stderr
    under_cov[:, i] = under_lower .< delta_true[2:end] .< under_upper

end

save(datadir("coverage/severe_case.jld2"), Dict(
    "correct_delta" => correct_delta,
    "under_delta" => under_delta,
    "correct_cov" => correct_cov,
    "under_cov" => under_cov,
    "delta_true" => delta_true
))

under_coef = inv(under_reg.omega) * under_reg.pi_mat
correct_coef = inv(correct_reg.omega) * correct_reg.pi_mat

perm_mat1 = perm_matrix(dimvals, under_rank)
perm_mat2 = perm_matrix(dimvals, true_rank)
perm_resp1 = (perm_mat1*cen_data)[:, 2:end]
perm_resp2 = (perm_mat2*cen_data)[:, 2:end]
pred = cen_data[:, 1:(end-1)]
resp1 = perm_resp1 .- mean(perm_resp1, dims=2)
resp2 = perm_resp2 .- mean(perm_resp2, dims=2)

res1 = resp1 - under_coef * pred
res2 = resp2 - correct_coef * pred

s1 = res1 * res1'
s2 = res2 * res2'

λvalues = 0.0:0.1:20.0
optimalλ1 = tuningselect(s1, 249, λvalues, verbose=false)
glasso1 = glasso(s1, obs, optimalλ1)
theta1 = glasso1.θ
heatmap(round.(theta1, digits=3), yflip=true)

optimalλ2 = tuningselect(s2, 249, λvalues, verbose=false)
glasso2 = glasso(s2, obs, 9)
theta2 = glasso2.θ
heatmap(round.(theta2, digits=3), yflip=true)

SSR1 = sum(res1 .^ 2)
SSR2 = sum(res2 .^ 2)




