
@testset "unpacking and repacking parameters" begin
    dimvals = [3, 4]
    ranks = [2, 1]
    p = 1
    theta = randn(92)
    delta, gamma, u3, u4, ll = unpack_params(theta, dimvals, ranks; p)

    @test size(u3) == (dimvals[1], ranks[1])
    @test size(u4) == (dimvals[2], ranks[2])
    @test size(delta) == (dimvals[1], dimvals[1] - ranks[1])
    @test size(gamma) == (dimvals[2], dimvals[2] - ranks[2])


    p = 2
    theta = randn(101)

    delta, gamma, u3, u4, ll = unpack_params(theta, dimvals, ranks; p)
    @test size(u3) == (p * dimvals[1], ranks[1])
    @test size(u4) == (p * dimvals[2], ranks[2])
end

@testset "unpacking and repacking parameters" begin
    dimvals = [3, 4]
    ranks = [2, 1]
    p = 2

    data = randn(12, 100)
    data = companion_data(data, p)
    perm_mat = kron(I(p), perm_matrix(dimvals, ranks))
    perm_resp = (perm_mat*data)[:, 2:end]
    pred = data[:, 1:(end-1)]
    resp = perm_resp .- mean(perm_resp, dims=2)
    iters = 100
    tol = 1e-05
    num_starts = 20
    num_selected = 10
    i = 1
end

