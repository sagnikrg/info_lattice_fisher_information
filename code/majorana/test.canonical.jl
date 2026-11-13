# checks canonical.jl against the free Majorana map at the Ising points, and the exact g=1 statements
using Test
include("canonical.jl"); include("majorana_floquet.jl")
@testset "canonical" begin
    L = 6; rng = MersenneTwister(3)
    g = 0.3; J = pi/4 .+ 0.2 .* rand(rng, L-1)
    U = floquet_kitaev(L; g = fill(g, L), Jz = J)
    Z1 = zdiag(1, L)
    Ct = [real(tr(U'^t * (Z1 .* U^t) * Diagonal(Z1))) / 2^L for t in 0:6]
    M = bond_map([J; 0.0], L; pbc = 0) * kick_map(g, L)
    Mt = [(M^t)[1, 1] for t in 0:6]
    M2 = kick_map(g, L) * bond_map([J; 0.0], L; pbc = 0)
    Mt2 = [(M2^t)[1, 1] for t in 0:6]
    @test isapprox(Ct, Mt; atol = 1e-10) || isapprox(Ct, Mt2; atol = 1e-10)
    # g = 1: every model anticommutes with Z_1 (pi pairing exact), incl. YY and cluster
    for (Jz, Jy, Kc) in ((J, zeros(L-1), zeros(L-2)), (J, 0.3 .* J, zeros(L-2)), (zeros(L-1), zeros(L-1), J[1:L-2]))
        U = floquet_kitaev(L; g = ones(L), Jz = Jz, Jy = Jy, Kc = Kc)
        Jy == zeros(L-1) && @test norm(U * Diagonal(Z1) + Diagonal(Z1) * U) < 1e-10 || norm(U' * Diagonal(Z1) * U + Diagonal(Z1)) < 1e-10
        om = angle.(eigvals(U))
        if Jy == zeros(L-1)
            @test median(antipodal(om)) < 1e-8
        else   # w != Delta: Z_1 is no longer exact, the dressed pi mode splits exponentially in L
            d = [median(antipodal(angle.(eigvals(floquet_kitaev(l; g = ones(l), Jz = fill(pi/4, l-1), Jy = fill(pi/32, l-1)))))) / 2^l for l in (6, 8, 10)]
            @test d[2] < d[1] / 10 && d[3] < d[1] / 10
        end
    end
    # psd: exact pairs -> odd line zero
    om = angle.(eigvals(floquet_kitaev(L; g = ones(L), Jz = J)))
    ev, od = psd_lines(om); @test od < 1e-12 * ev
end
