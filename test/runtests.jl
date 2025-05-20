using DrWatson
@quickactivate :PseudoStructuralComovements

println("Starting tests")
ti = time()

@testset "PseudoStructuralComovements tests" begin
    Random.seed!(20250421)

    dimvals = [3, 4]
    ranks = [2, 3]

    coef = generate_rrmar_coef(dimvals, ranks)
    u3_true = coef.u3
    u4_true = coef.u4

    delta_true = coef.delta
    gamma_true = coef.gamma

    delta_star = delta_true[2:end]
    gamma_star = gamma_true[2:end]

    omega_true = omega_from_both(delta_star, gamma_star, dimvals, ranks)
    pi_true = pi_from_both(u3_true, u4_true, dimvals, ranks)
    perm_mat = both_perm_mat(dimvals, ranks)

    constructed_coef = inv(omega_true) * pi_true
    coef.A
    kron(u2_true, u1_true) * kron(u4_true, u3_true)'
    norm(coef.A - constructed_coef)

end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60, digits=3), " minutes")

using LinearAlgebra, Statistics, TensorToolbox

N1_r1 = dimvals[1] - ranks[1]
N2_r2 = dimvals[2] - ranks[2]
r = prod(ranks)
pdims = prod(dimvals)

coef = generate_rrmar_coef(dimvals, ranks)

tensor_coef =
    matten(coef.A, [1, 2], [3, 4], [dimvals[1], dimvals[2], dimvals[1], dimvals[2]])
flat1 = tenmat(tensor_coef, row=1)
flat2 = tenmat(tensor_coef, row=2)
flat3 = tenmat(tensor_coef, row=3)
flat4 = tenmat(tensor_coef, row=4)
u1 = svd(flat1).U[:, 1:ranks[1]]
u2 = svd(flat2).U[:, 1:ranks[2]]

delta = nullspace(u1')
rotate_u!(delta)
delta_star = delta[(N1_r1+1):end, :]
gamma = nullspace(u2')
rotate_u!(gamma)
gamma_star = gamma[(N2_r2+1):end, :]

perm_mat = both_perm_mat(dimvals, ranks)
omega = omega_from_both(delta_star, gamma_star, dimvals, ranks)
pi_mat = omega * perm_mat * coef.A
pi_star = pi_mat[(pdims-r+1):end, :]
u4_est, u3_est = nearest_kron(pi_star', size(u2), size(u1))

s = u3_est[1, 1]

u3 = u3_est / s
u4 = u4_est * s





