# Sanity checks for gh_probe.jl, L=4..6, seconds.
using Test
include("gh_probe.jl")

@testset "gh_probe" begin
    L = 6; x0 = 3
    M = gh_initial_state(L, x0, rand(0:1, L))
    @test norm(M) ≈ 1
    # Bell pair: S(R)=1, I(R:x0)=2, I(R:j)=0 elsewhere, <ZZ>=<XX>=1, <YY>=-1
    @test entropy_bits(reduced_dm(M, true, Int[], L)) ≈ 1 atol=1e-10
    @test mutual_info_R(M, [x0], L) ≈ 2 atol=1e-10
    @test all(abs(mutual_info_R(M, [j], L)) < 1e-10 for j in 1:L if j != x0)
    @test R_site_correlator(M, PZ, x0, L) ≈ 1
    @test R_site_correlator(M, PX, x0, L) ≈ 1
    @test R_site_correlator(M, PY, x0, L) ≈ -1
    # the whole chain always has I(R:chain)=2, S(R)=1 under the circuit
    r = gh_probe_run(L, 0.3, 0.5, 20; x0 = x0)
    @test all(abs.(r.SR .- 1) .< 1e-8)
    @test all(abs.(r.Iblock[:, L] .- 2) .< 1e-8)
    # g=0, theta=0: U diagonal in Z up to the small Z/theta noise, so the ZZ correlator stays near 1
    r0 = gh_probe_run(L, 0.0, 0.0, 10; x0 = x0)
    @test minimum(r0.ZZ) > 0.9
    # g=1 perfect kick flips Z every period: ZZ alternates in sign
    r1 = gh_probe_run(L, 0.0, 1.0, 10; x0 = x0)
    @test all(r1.ZZ[2:2:end] .< -0.5) && all(r1.ZZ[3:2:end] .> 0.5)
end
