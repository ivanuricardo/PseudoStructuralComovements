using DrWatson
@quickactivate :PseudoStructuralComovements

matdata = load(datadir("./country_indicators/globaldata.jld2"), "matdata");
vecdata = tenmat(matdata, row=[1, 2])
cen_data = vecdata .- mean(vecdata, dims=2)

dimvals = collect(size(matdata[:, :, 1]))
ranks = [2,1]

#=res = comovement_reg(cen_data, dimvals, ranks)=#

#=save(datadir("country_indicators/results.jld2"), Dict(=#
#=    "res" => res,=#
#=))=#

results_dict = load(datadir("country_indicators/results.jld2"))
res = results_dict["res"]





res.delta_est
res.delta_stderr
res.gamma_est
res.gamma_stderr

delta_star = res.delta_est[2:end]
gamma_star = res.gamma_est[end, :]'
joint_kron = kron(gamma_star, delta_star)'

inv_hessian = inv(res.hess_est)[1:5, 1:5]
sqrt.(diag(inv_hessian))
# We need the covariances

c11 = inv_hessian[3,1]
c12 = inv_hessian[4,1]
c13 = inv_hessian[5,1]
c21 = inv_hessian[3,2]
c22 = inv_hessian[4,2]
c23 = inv_hessian[5,2]

d1_g1_std = sqrt(delta_star[1]^2 * res.gamma_stderr[1]^2 + gamma_star[1]^2 * res.delta_stderr[1]^2 + 2 * delta_star[1] * gamma_star[1] * c11)
d1_g2_std = sqrt(delta_star[1]^2 * res.gamma_stderr[2]^2 + gamma_star[2]^2 * res.delta_stderr[1]^2 + 2 * delta_star[1] * gamma_star[2] * c12)
d1_g3_std = sqrt(delta_star[1]^2 * res.gamma_stderr[3]^2 + gamma_star[3]^2 * res.delta_stderr[1]^2 + 2 * delta_star[1] * gamma_star[3] * c13)

d2_g1_std = sqrt(delta_star[2]^2 * res.gamma_stderr[1]^2 + gamma_star[1]^2 * res.delta_stderr[2]^2 + 2 * delta_star[2] * gamma_star[1] * c21)
d2_g2_std = sqrt(delta_star[2]^2 * res.gamma_stderr[2]^2 + gamma_star[2]^2 * res.delta_stderr[2]^2 + 2 * delta_star[2] * gamma_star[2] * c22)
d2_g3_std = sqrt(delta_star[2]^2 * res.gamma_stderr[3]^2 + gamma_star[3]^2 * res.delta_stderr[2]^2 + 2 * delta_star[2] * gamma_star[3] * c23)
