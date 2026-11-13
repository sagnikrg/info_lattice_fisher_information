using Test, Random
include("ed_circuit.jl"); include("canonical.jl"); include("gh_edge.jl")
@testset "gh_edge" begin
    L = 6; rng = MersenneTwister(1)
    # perfect point g=1, theta=0, any fields: Z_1 exact pi mode -> Api_ZZ = 1 at the edge
    U = circuit_ctrl(L, 0.0, 1.0, rand(rng, L) * 2pi, rand(rng, L - 1) * pi)
    r = gh_edge_run(U, L, 1, rand(rng, 0:1, L); t1 = 10, t2 = 30)
    @test r.Api_ZZ ≈ 1 atol = 1e-10
    @test r.A0_ZZ ≈ 0 atol = 0.06          # odd window length leaves 1/n
    # g=0: Z conserved everywhere -> A0_ZZ = 1 and I(R:x0) = 2 at any site... ZZ sign constant
    U0 = circuit_ctrl(L, 0.0, 0.0, rand(rng, L) * 2pi, rand(rng, L - 1) * pi)
    r0 = gh_edge_run(U0, L, 3, rand(rng, 0:1, L); t1 = 10, t2 = 30)
    @test r0.A0_ZZ ≈ 1 atol = 1e-10
    # ZZ correlator equals the infinite temperature like autocorrelator for a single config: |.|<=1
    @test abs(r.Api_XX) <= 1 + 1e-12
end
