using DrWatson, Test
@quickactivate "PseudoStructuralComovements"

println("Starting tests")
ti = time()

@testset "PseudoStructuralComovements tests" begin
    Random.seed!(20250421)

    dimvals = [3, 4]
    true_ranks = [2, 3]
    under_rank = [1, 3]
    over_rank = [3, 3]

    sims = 50
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
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60, digits=3), " minutes")
