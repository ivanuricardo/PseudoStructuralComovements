module PseudoStructuralComovements

using Reexport
@reexport using LinearAlgebra, Statistics, Random, JLD2, StatsPlots
@reexport using ProgressMeter

using TensorToolbox
using Optim
using SparseArrays
using ForwardDiff
using ReverseDiff
using Distributions
using NLSolversBase
using ProximalOperators

include("blikelihood.jl")
export b_unpack_params
export b_pack_params
export init_both
export rand_init
export comovement_reg

include("./bothfuncs.jl")
export both_blockvec
export both_perm_mat
export b_top_left
export b_right1
export b_right2
export b_right3
export omega_from_both
export pi_from_both
export both_loglike
export b_secondrow
export b_thirdrow
export b_fourthrow

include("./helper.jl")
export insertk!
export removek!
export cov_to_ll
export make_companion
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

include("simtools.jl")
export isstable
export generate_rrmar_coef
export simulate_rrmar_data
export sim_stats

end
