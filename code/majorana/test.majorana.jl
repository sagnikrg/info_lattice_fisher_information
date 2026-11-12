using Test, Random
include("ed_circuit.jl")
include("majorana_floquet.jl")
Random.seed!(1)
@testset "JW map = ED at theta=0, h=0, no noise" begin
    for (L, g) in ((6, 0.9), (6, 0.3), (8, 0.97))
        J = rand(L - 1) * pi
        U = circuit_ctrl(L, 0.0, g, zeros(L), J)
        ed = edge_autocorr_ed(U, 1, 20, L)
        M = floquet_majorana(g, J)
        ff = [(M^t)[1, 1] for t in 0:20]
        @test maximum(abs.(ed .- ff)) < 1e-10
        # last site: Z_L = (prod X) g_{2L-1}, its autocorrelator is not a single Majorana one; skip
    end
end
@testset "pfaffian sign" begin
    for n in (2, 4, 6, 8)
        A = randn(n, n); A = A - A'
        @test pfaffian_sign(A)^2 == 1
        @test abs(det(A)) > 0
    end
    A = [0 1.0; -1 0]; @test pfaffian_sign(A) == 1
    # Pf of block diag = product
    B = zeros(4, 4); B[1, 2] = 2; B[2, 1] = -2; B[3, 4] = -3; B[4, 3] = 3
    @test pfaffian_sign(B) == -1
end
