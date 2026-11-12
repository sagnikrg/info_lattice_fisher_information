using Random, Printf, Statistics
include("majorana_floquet.jl")
# Pf brute check n=4
A=randn(4,4); A=A-A'; pf=A[1,2]*A[3,4]-A[1,3]*A[2,4]+A[1,4]*A[2,3]; @assert sign(pf)==pfaffian_sign(A)
Random.seed!(2)
L = 100
println("clean J=pi/4 (bond rotation pi/2) and J=0.3:")
for Jc in (pi/4, 0.3), g in 0.0:0.1:1.0
    g == 0 && continue
    Jr = fill(Jc, L); Jo = fill(Jc, L-1)
    nu = floquet_invariants(g, Jr); em = edge_modes(floquet_majorana(g, Jo); tol=1e-6)
    @printf("J=%.3f g=%.2f  nu0=%d nupi=%d | open chain modes at 0: %d  at pi: %d\n", Jc, g, nu.nu0, nu.nupi, em.n0, em.npi)
end
println("random J in [0,pi], 50 realisations, fractions:")
for g in 0.05:0.1:0.95
    r = [begin J = rand(L)*pi; nu = floquet_invariants(g, J); em = edge_modes(floquet_majorana(g, J[1:L-1]); tol=1e-6); (nu.nu0, nu.nupi, em.n0>0, em.npi>0) end for _ in 1:50]
    @printf("g=%.2f  <nu0>=%.2f <nupi>=%.2f | P(0 mode)=%.2f P(pi mode)=%.2f\n", g, mean(first.(r)), mean(getindex.(r,2)), mean(getindex.(r,3)), mean(getindex.(r,4)))
end
