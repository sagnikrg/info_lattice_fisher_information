# Free Majorana chain (theta=0, h=0, J uniform in [0,pi]): chiral invariants on the ring
# and 0/pi edge mode splittings on the open chain, g scan, L = 32, 64, 128.
using Random, HDF5, Statistics, LinearAlgebra
BLAS.set_num_threads(1)
include("majorana_floquet.jl")
out = joinpath(@__DIR__, "..", "..", "data", "majorana"); mkpath(out)
fname = joinpath(out, "free_scan.hdf5")
gs = collect(0.0:0.025:1.0); NR = 200
Random.seed!(20261112)
for L in (32, 64, 128)
    nu0 = fill(NaN, length(gs), NR); nupi = similar(nu0); d0 = similar(nu0); dpi = similar(nu0)
    for (ig, g) in enumerate(gs), r in 1:NR
        J = rand(L) * pi
        c = try chiral_invariants(g, J) catch; (nu0 = NaN, nupi = NaN) end
        nu0[ig, r] = c.nu0; nupi[ig, r] = c.nupi
        q = abs.(angle.(eigvals(floquet_majorana(g, J[1:L-1]))))   # |quasienergy| in [0, pi]
        d0[ig, r] = minimum(q); dpi[ig, r] = pi - maximum(q)          # splittings of the 0 and pi pairs
    end
    h5open(fname, "cw") do f
        for (k, v) in (("nu0", nu0), ("nupi", nupi), ("split0", d0), ("splitpi", dpi))
            haskey(f, "L$L/$k") && delete_object(f, "L$L/$k"); f["L$L/$k"] = v
        end
        haskey(f, "g") || (f["g"] = gs)
    end
    println("L=$L done"); flush(stdout)
end
