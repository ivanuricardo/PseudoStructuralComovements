module PseudoStructuralComovements

const R_LOCK = ReentrantLock()
export R_LOCK

using Reexport
@reexport using LinearAlgebra, Statistics, Random, JLD2, StatsPlots
@reexport using ProgressMeter, TensorToolbox, Makie, CairoMakie

using Optim
using RCall
using SparseArrays
using Distributions
using NLSolversBase

include("likelihood.jl")
export unpack_params
export pack_params
export init_alg
export loglike
export rand_init
export comovement_reg
export rrmar

include("./likelihood_helpers.jl")
export blockvec
export perm_matrix
export top_left
export right1
export right2
export right3
export create_omega
export create_pi
export create_second_row
export create_third_row
export create_fourth_row

include("./helper.jl")
export insertk!
export check_rank
export removek!
export cov_to_ll
export make_companion
export companion_data
export rotate_u!
export ll_to_cov
export vec_to_ll
export vecb
export vech
export check_conf
export var_coef
export ols_coef
export nearest_kron
export nearest_posdef

include("informationcrit.jl")
export system_parameters
export aic
export bic
export hqc
export rank_selection

include("simulation_helpers.jl")
export isstable
export generate_rrmar_coef
export simulate_rrmar_data
export sim_stats

function __init__()
    r_path1 = joinpath(@__DIR__, "r_helpers.R")  # Path to R file
    r_path2 = joinpath(@__DIR__, "tenAR.R")  # Path to R file
    r_path3 = joinpath(@__DIR__, "tenFM.R")  # Path to R file
    R"source($r_path1)"  # Source the file in R
    R"source($r_path2)"  # Source the file in R
    R"source($r_path3)"  # Source the file in R
end

end
