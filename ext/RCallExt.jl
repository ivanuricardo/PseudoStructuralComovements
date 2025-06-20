module RCallExt

# only runs if the user has added RCall to their environment
using PseudoStructuralComovements
using RCall

function __init__()
    if Sys.which("R") !== nothing
        R"""
        if (!requireNamespace("tensorTS", quietly=TRUE))
            install.packages("tensorTS", repos="https://cloud.r-project.org")
        """
    else
        @warn "R binary not found—skipping tensorTS installation"
    end
end

# now extend your main package:
PseudoStructuralComovements.run_r_model(args...) = begin
    #= call into R using RCall=#
    #=RCall.reval("some_R_function(...)")=#
end

end # module
