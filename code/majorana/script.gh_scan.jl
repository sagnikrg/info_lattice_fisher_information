# Reference qubit probe as an order parameter: g scans of our circuit (fields h_j uniform in [0,2pi), J_j uniform in [0,pi), noise off)
# at theta = 0, 0.1, 0.2, 0.3, 0.4, R on the edge qubit and on the bulk qubit L/2, L = 8, 10.
# usage: julia -t 6 script.gh_scan.jl  -> data/majorana/gh_scan.hdf5 (per realisation late time signals), log_gh_scan.txt
using HDF5, Random, Statistics, LinearAlgebra
BLAS.set_num_threads(1)
include("ed_circuit.jl"); include("gh_edge.jl")
OUT = joinpath(@__DIR__, "..", "..", "data", "majorana"); fname = joinpath(OUT, "gh_scan.hdf5")
gs = collect(0.0:0.05:1.0); ths = (0.0, 0.1, 0.2, 0.3, 0.4); NR = Dict(8 => 50, 10 => 30)
keysig = (:I_R, :Api_ZZ, :A0_ZZ)
logf(s) = (open(joinpath(OUT, "log_gh_scan.txt"), "a") do io; println(io, s); end; println(s))
for L in (8, 10), th in ths
    grp = "L$L/theta$th"
    isfile(fname) && h5open(f -> haskey(f, grp), fname, "r") && continue
    res = Dict((x, k) => zeros(length(gs), NR[L]) for x in (:edge, :bulk), k in keysig)
    t = @elapsed Threads.@threads for idx in 1:length(gs) * NR[L]
        ig, r = divrem(idx - 1, NR[L]) .+ 1
        rng = Xoshiro(hash((L, th, ig, r, :gh_scan)))
        U = circuit_ctrl(L, th, gs[ig], rand(rng, L) * 2pi, rand(rng, L - 1) * pi); cfg = rand(rng, 0:1, L)
        for (x, x0) in ((:edge, 1), (:bulk, L ÷ 2))
            s = gh_edge_run(U, L, x0, cfg)
            for k in keysig; res[(x, k)][ig, r] = getfield(s, k); end
        end
    end
    h5open(fname, "cw") do f
        for ((x, k), v) in res; f["$grp/$x/$k"] = v; end
        haskey(f, "g") || (f["g"] = gs)
    end
    m(x, k) = round.(vec(mean(res[(x, k)], dims = 2))[1:4:end], digits = 2)
    logf("$grp  edge pi $(m(:edge,:Api_ZZ))  bulk pi $(m(:bulk,:Api_ZZ))  bulk 0 $(m(:bulk,:A0_ZZ))  $(round(t, digits=1))s")
end
