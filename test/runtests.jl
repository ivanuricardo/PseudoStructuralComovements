using DrWatson, Test
@quickactivate :PseudoStructuralComovements

println("Starting tests")
ti = time()

@testset "PseudoStructuralComovements basic test" begin
    @test 1 == 1
end
include("./test-likelihood.jl")

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60, digits=3), " minutes")
