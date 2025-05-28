
@testset "unpacking and repacking parameters" begin
    dimvals = [3, 4]
    ranks = [2, 1]
    p = 1
    theta = randn(92)
    delta, gamma, u3, u4, ll = b_unpack_params(theta, dimvals, ranks; p)

    @test size(u3) == (dimvals[1], ranks[1])
    @test size(u4) == (dimvals[2], ranks[2])
    @test size(delta) == (dimvals[1], dimvals[1] - ranks[1])
    @test size(gamma) == (dimvals[2], dimvals[2] - ranks[2])


    theta_repacked = b_pack_params(delta, gamma, u3, u4, ll; p)

    p = 2
    theta = randn(101)

    delta, gamma, u3, u4, ll = b_unpack_params(theta, dimvals, ranks; p)
    @test size(u3) == (p * dimvals[1], ranks[1])
    @test size(u4) == (p * dimvals[2], ranks[2])
end
