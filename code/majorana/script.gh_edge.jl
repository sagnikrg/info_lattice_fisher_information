# Gullans-Huse probe on the pi mode: edge (x0=1) vs bulk (x0=L/2) reference qubit, our circuit and the canonical chains.
# usage: julia -t 1 script.gh_edge.jl   -> data/majorana/gh_edge.hdf5 (per realisation late time signals, small), log_gh_edge.txt
using HDF5, Random, Statistics
include("ed_circuit.jl"); include("canonical.jl"); include("gh_edge.jl")
OUT = joinpath(@__DIR__, "..", "..", "data", "majorana"); fname = joinpath(OUT, "gh_edge.hdf5")
NR = Dict(8 => 50, 10 => 30)
ours = [(:dtc, 0.0, 0.95), (:dtc_theta, 0.2, 0.95), (:pm, 0.0, 0.05), (:scrambling, 0.35, 0.5)]
canon = [:SG0, :PM, :SGpi, :SGpi_strong, :PM0pi, :cluster_pi]
keysig = (:I_R, :Api_ZZ, :A0_ZZ, :Api_XX, :A0_XX, :Api_YY, :A0_YY)
logf(s) = (open(joinpath(OUT, "log_gh_edge.txt"), "a") do io; println(io, s); end; println(s))
for L in (8, 10)
    pts = [[(nm, rng -> circuit_ctrl(L, th, g, rand(rng, L) * 2pi, rand(rng, L - 1) * pi)) for (nm, th, g) in ours];
           [(nm, rng -> realisation(POINTS[findfirst(p -> p[1] == nm, POINTS)], L, rng)) for nm in canon]]
    for (ip, (nm, mk)) in enumerate(pts)
        grp = "L$L/$nm"
        isfile(fname) && h5open(f -> haskey(f, grp), fname, "r") && continue
        rng = MersenneTwister(7000L + ip)
        res = Dict((x, k) => zeros(NR[L]) for x in (:edge, :bulk), k in keysig)
        t = @elapsed for r in 1:NR[L]
            U = mk(rng); cfg = rand(rng, 0:1, L)
            for (x, x0) in ((:edge, 1), (:bulk, L ÷ 2))
                s = gh_edge_run(U, L, x0, cfg)
                for k in keysig; res[(x, k)][r] = getfield(s, k); end
            end
        end
        h5open(fname, "cw") do f
            for ((x, k), v) in res; f["$grp/$x/$k"] = v; end
        end
        m(x, k) = round(mean(res[(x, k)]), digits = 3)
        logf("L$L $(rpad(string(nm), 12)) edge I=$(m(:edge,:I_R)) piZZ=$(m(:edge,:Api_ZZ)) 0ZZ=$(m(:edge,:A0_ZZ)) piXX=$(m(:edge,:Api_XX)) | bulk I=$(m(:bulk,:I_R)) piZZ=$(m(:bulk,:Api_ZZ)) 0ZZ=$(m(:bulk,:A0_ZZ)) 0XX=$(m(:bulk,:A0_XX))  $(round(t, digits=1))s")
    end
end
