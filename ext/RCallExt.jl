module RCallExt

using PseudoStructuralComovements
using RCall

#=function __init__()=#
#=    if Sys.which("R") !== nothing=#
#=        R"""=#
#=        if (!requireNamespace("tensorTS", quietly=TRUE)) {=#
#=            install.packages("tensorTS", repos="https://cloud.r-project.org")}=#
#=        library(tensorTS)=#
#=        RCall.reval(string("source(\"", joinpath(@__DIR__, "rfuncs.R"), "\")"))=#
#=        """=#
#=    else=#
#=        @warn "R binary not found—skipping tensorTS installation"=#
#=    end=#
#=end=#

function run_r_regression(x::Vector{Float64}, y::Vector{Float64})
    @rput x y
    R"""
    model <- lm(y ~ x)
    coef(model)
    """
end

export run_r_regression

end
