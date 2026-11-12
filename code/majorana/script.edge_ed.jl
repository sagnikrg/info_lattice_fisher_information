# Infinite temperature edge autocorrelator C(t) = <Z_1(t) Z_1>, ED of the thesis circuit with
# noise off. Late time signals: pi mode  A_pi = mean_t (-1)^t C(t),  0 mode A_0 = mean_t C(t),
# t in 200..400, averaged over realisations. Cases: h = 0 or random in [0,2pi]; theta = 0, 0.1, 0.2.
using Random, HDF5, Statistics, LinearAlgebra
BLAS.set_num_threads(1)
include("ed_circuit.jl")
out = joinpath(@__DIR__, "..", "..", "data", "majorana"); mkpath(out)
fname = joinpath(out, "edge_ed.hdf5")
gs = [0.0, 0.05, 0.1, 0.2, 0.3, 0.7, 0.8, 0.9, 0.95, 1.0]; T = 400; tl = 201:401
Random.seed!(7)
for L in (8, 10), hon in (false, true), theta in (0.0, 0.1, 0.2, 0.3)
    NR = L == 10 ? 16 : 100
    Api = zeros(length(gs), NR); A0 = similar(Api)
    for (ig, g) in enumerate(gs), r in 1:NR
        J = rand(L - 1) * pi; h = hon ? rand(L) * 2pi : zeros(L)
        C = edge_autocorr_ed(circuit_ctrl(L, theta, g, h, J), 1, T, L)
        Api[ig, r] = mean(((-1) .^ (tl .- 1)) .* C[tl]); A0[ig, r] = mean(C[tl])
    end
    grp = "L$L/h$(Int(hon))/theta$theta"
    h5open(fname, "cw") do f
        haskey(f, grp) && delete_object(f, grp)
        f["$grp/Api"] = Api; f["$grp/A0"] = A0; haskey(f, "g") || (f["g"] = gs)
    end
    println(grp, "  A_pi: ", round.(vec(mean(Api, dims=2)), digits=3), "  A_0: ", round.(vec(mean(A0, dims=2)), digits=3)); flush(stdout)
end
