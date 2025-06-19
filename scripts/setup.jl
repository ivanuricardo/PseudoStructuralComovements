using DrWatson
@quickactivate :PseudoStructuralComovements

if Sys.which("R") !== nothing
    R"""
    if (!requireNamespace("tensorTS", quietly = TRUE)) {
      install.packages("tensorTS", repos = "https://cloud.r-project.org")
    }
    """
else
    @warn "R is not installed or not found in PATH. Skipping R-related operations."
end
