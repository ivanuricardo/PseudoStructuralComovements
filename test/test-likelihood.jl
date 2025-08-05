
@testset "unpacking and repacking parameters" begin
    dimvals = [3, 4]
    ranks = [2, 1]
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]
    p = 1
    theta = randn(29)
    delta, gamma, u3, u4, ll1, ll2 = unpack_params(theta, dimvals, ranks; p)

    @test size(u3) == (dimvals[1], ranks[1])
    @test size(u4) == (dimvals[2], ranks[2])
    @test size(delta) == (dimvals[1], dimvals[1] - ranks[1])
    @test size(gamma) == (dimvals[2], dimvals[2] - ranks[2])

    # Checking if repack gives me back theta
    delta_star = delta[(N1_r1+1):end, :]
    gamma_star = gamma[(N2_r2+1):end, :]
    theta_reconstructed = pack_params(delta_star, gamma_star, u3, u4, ll1, ll2)
    @test norm(theta_reconstructed - theta) ≈ 0


    # Testing with 2 lags
    p = 2
    theta = randn(38)

    delta, gamma, u3, u4, ll1, ll2 = unpack_params(theta, dimvals, ranks; p)
    @test size(u3) == (p * dimvals[1], ranks[1])
    @test size(u4) == (p * dimvals[2], ranks[2])

    # repacking again should give me the same
    delta_star = delta[(N1_r1+1):end, :]
    gamma_star = gamma[(N2_r2+1):end, :]
    theta_reconstructed = pack_params(delta_star, gamma_star, u3, u4, ll1, ll2; p)
    @test norm(theta_reconstructed - theta) ≈ 0

    # Testing with 3 lags
    p = 3
    theta = randn(47)

    delta, gamma, u3, u4, ll1, ll2 = unpack_params(theta, dimvals, ranks; p)
    @test size(u3) == (p * dimvals[1], ranks[1])
    @test size(u4) == (p * dimvals[2], ranks[2])

    # repacking again should give me the same
    delta_star = delta[(N1_r1+1):end, :]
    gamma_star = gamma[(N2_r2+1):end, :]
    theta_reconstructed = pack_params(delta_star, gamma_star, u3, u4, ll1, ll2; p)
    @test norm(theta_reconstructed - theta) ≈ 0
end

@testset "check permutation matrix" begin
    Y = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]
    dimvals = [4,4]
    ranks = [1,3]
    # Y11 should be dimension 3x1, Y12 should be of dimension 3x3
    # Y21 should be of dimension 1x1, Y22 should be of dimension 1x3
    bvec_y = blockvec(Y, dimvals, ranks)
    Y11 = bvec_y[1:3] 
    @test isapprox(Y11, [1,5,9])
    Y12 = bvec_y[4]
    @test isapprox(Y12, 13)
    Y21 = bvec_y[5:13]
    @test isapprox(Y21, [2,6,10,3,7,11,4,8,12])
    Y22 = bvec_y[14:end]
    @test isapprox(Y22, [14,15,16])

    perm_mat = perm_matrix(dimvals, ranks)
    permed_y = perm_mat * vec(Y)
    @test isapprox(permed_y, bvec_y)
    # P y is our output

end

@testset "omega and pi matrices create a" begin
    dimvals = [3,4]
    ranks = [2,1]
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]
    p = 1
    theta = randn(29)
    delta, gamma, u3, u4, ll1, ll2 = unpack_params(theta, dimvals, ranks; p)
    delta_star = delta[(N1_r1+1):end, :]
    gamma_star = gamma[(N2_r2+1):end, :]
    omega = create_omega(delta_star, gamma_star, dimvals, ranks)
    pi = create_pi(u3, u4, dimvals, ranks; p=1)
    u1 = vcat(-delta_star', I)
    u2 = vcat(-gamma_star', I)
    perm_mat = perm_matrix(dimvals, ranks)

    A = kron(u2, u1) * kron(u4, u3)'
    omega_pi = perm_mat' * inv(omega) * pi
    # Should be P' inv(Ω^*) Π = A. In other words, P' inv(Ω P) Π = A
    # This is the form that should be used in likelihood, all cases
    @test isapprox(norm(A), norm(omega_pi))
end
