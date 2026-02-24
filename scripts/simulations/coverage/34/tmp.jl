using DrWatson
@quickactivate :PseudoStructuralComovements
Random.seed!(20250723)

dimvals = [2, 3]
n = prod(dimvals)
true_rank = [1, 1]

burnin = 200
obs = 300 + burnin

coef = generate_rrmar_coef(dimvals, true_rank)
delta_true = coef.delta
delta_star = delta_true[2:end]
gamma_true = coef.gamma
gamma_star = gamma_true[3, :]'
u3_true = coef.u3
u4_true = coef.u4

u1_true = nullspace(delta_true') * inv(nullspace(delta_true')[2])
u2_true = nullspace(gamma_true') * inv(nullspace(gamma_true')[3])

norm(coef.A - kron(u2_true, u1_true) * kron(u4_true, u3_true)')

omega = create_omega(delta_star, gamma_star, dimvals, true_rank)
pi_mat = create_pi(u3_true, u4_true, dimvals, true_rank)

inv(omega) * pi_mat

data = simulate_rrmar_data(dimvals, true_rank, obs; A=coef, burnin, matrix_err=true)
correct_reg = comovement_reg(data.data, dimvals, true_rank; iters=1000)
inv(correct_reg.omega) * correct_reg.pi_mat

# Permute the rows of the identity to obtain a switching of the rows
P1 = I(dimvals[1])[randperm(dimvals[1]), :]
P2 = I(dimvals[2])[randperm(dimvals[2]), :]
P = kron(P2, P1)

perm_reg = comovement_reg(P * data.data, dimvals, true_rank; iters = 1000)

perm_est = P * inv(perm_reg.omega) * perm_reg.pi_mat * P'
est = inv(correct_reg.omega) * correct_reg.pi_mat
norm(perm_est - est)


norm(perm_est - est)

perm_data = P * data.data
perm_resp = perm_data[:, 2:end] .- mean(perm_data[:, 2:end], dims = 2)
perm_pred = perm_data[:, 1:end-1] .- mean(perm_data[:, 1:end-1], dims = 2)
perm_coef = copy(P * est * P')
perm_res = perm_resp - perm_coef * perm_pred
perm_cov = perm_res * perm_res' / 299
logdet(perm_cov)


resp = data.data[:, 2:end] .- mean(data.data[:, 2:end], dims = 2)
pred = data.data[:, 1:end-1] .- mean(data.data[:, 1:end-1], dims = 2)
coef = inv(omega) * pi_mat
res = resp - coef * pred
cov_mat = res * res' / 299
logdet(cov_mat)


