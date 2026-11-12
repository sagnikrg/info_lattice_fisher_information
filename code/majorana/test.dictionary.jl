using Test, LinearAlgebra
include("ed_circuit.jl")   # X, Y, Z, site_op
L = 4
Xs(k) = site_op(ComplexF64.(X), k, L); Ys(k) = site_op(ComplexF64.(Y), k, L); Zs(k) = site_op(ComplexF64.(Z), k, L)
string_(j) = j == 1 ? Matrix{ComplexF64}(I, 2^L, 2^L) : prod(Xs(k) for k in 1:j-1)
gam(a) = isodd(a) ? string_((a + 1) ÷ 2) * Zs((a + 1) ÷ 2) : string_(a ÷ 2) * Ys(a ÷ 2)
@testset "JW dictionary" begin
    for a in 1:2L, b in 1:2L
        @test gam(a) * gam(b) + gam(b) * gam(a) ≈ 2 * (a == b) * I
    end
    for j in 1:L
        @test Xs(j) ≈ im * gam(2j - 1) * gam(2j)
        @test Zs(j) ≈ im^(j - 1) * prod(gam(a) for a in 1:2j-1)
    end
    for j in 1:L-1
        a, b, c, d = 2j - 1, 2j, 2j + 1, 2j + 2
        @test Zs(j) * Zs(j + 1) ≈ im * gam(b) * gam(c)
        @test Xs(j) * Xs(j + 1) ≈ -gam(a) * gam(b) * gam(c) * gam(d)
        @test Ys(j) * Ys(j + 1) ≈ -im * gam(a) * gam(d)
    end
    @test prod(Xs(j) for j in 1:L) ≈ im^L * prod(gam(a) for a in 1:2L)
end
