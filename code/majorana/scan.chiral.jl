using Random, Printf, Statistics
include("majorana_floquet.jl")
Random.seed!(3); L = 120
# chirality check
let J = rand(L) * pi, Kh = kick_map(0.4, L), Dh = bond_map(J ./ 2, L), G = Diagonal([isodd(i) ? 1.0 : -1.0 for i in 1:2L])
    M1 = Kh * Dh * Kh; println("chiral check |G M G - M^T| = ", norm(G * M1 * G - M1'))
end
for Jc in (pi/4, 0.3), g in 0.1:0.1:1.0
    c = chiral_invariants(g, fill(Jc, L)); em = edge_modes(floquet_majorana(g, fill(Jc, L-1)); tol=1e-6)
    @printf("J=%.3f g=%.1f nu1=%6.3f nu2=%6.3f -> nu0=%6.3f nupi=%6.3f | modes 0:%d pi:%d\n", Jc, g, c.nu1, c.nu2, c.nu0, c.nupi, em.n0, em.npi)
end
println("random J in [0,pi]:")
for g in 0.05:0.1:0.95
    r = [begin J = rand(L)*pi; c = chiral_invariants(g, J); em = edge_modes(floquet_majorana(g, J[1:L-1]); tol=1e-6); (c.nu0, c.nupi, em.n0>0, em.npi>0) end for _ in 1:40]
    @printf("g=%.2f <nu0>=%6.3f <nupi>=%6.3f | P(0)=%.2f P(pi)=%.2f\n", g, mean(first.(r)), mean(getindex.(r,2)), mean(getindex.(r,3)), mean(getindex.(r,4)))
end
