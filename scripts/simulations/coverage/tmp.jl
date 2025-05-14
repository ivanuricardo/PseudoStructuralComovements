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

data = simulate_rrmar_data(dimvals, true_ranks, obs; A=coef, burnin)
cen_data = data.data .- mean(data.data, dims=2)

correct_reg = comovement_reg(cen_data, dimvals, true_ranks, iters=100)
over_reg = comovement_reg(cen_data, dimvals, over_rank, iters=100)
under_reg = comovement_reg(cen_data, dimvals, under_rank, iters=100)


