#####################################################################
# Small test of the Gullans-Huse reference qubit probe (K.A.R.E.N., Ch5)
# usage: julia -t 1 script.gh_probe_small.jl L NREAL T
# The three example points of the thesis (PM theta=0 g=0.05, scrambling
# theta=0.35 g=0.5, DTC theta=0 g=0.95) and a coarse g scan at theta=0.
# Output: ../../data/gh_probe/gh_probe_small_L<L>.hdf5, per point the
# realisation arrays (SR, Isite, Iblock, ZZ, XX, YY) and the disorder.
#####################################################################

using HDF5, Dates, Random
include("gh_probe.jl")

L     = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 8
NREAL = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 50
T     = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 100
x0    = L ÷ 2

points = [(0.0, 0.05), (0.35, 0.5), (0.0, 0.95)]                       # examples
append!(points, [(0.0, g) for g in 0.1:0.1:0.9])                        # g scan, theta=0

outdir = joinpath(@__DIR__, "..", "..", "data", "gh_probe")
mkpath(outdir)
fname = joinpath(outdir, "gh_probe_small_L$(L).hdf5")

h5open(fname, "cw") do f
    a = HDF5.attributes(f)
    a["model"] = "circuit_dtc(L, theta, 1-g, h, J) from .header circuit/time_crystals.jl; g=1 perfect pi kick"
    a["probe"] = "reference R in a Bell pair with site x0, other sites a random Z product state; no measurements; S(R)=1 by construction"
    a["x0"] = x0; a["T"] = T; a["NREAL"] = NREAL
    a["code_probe"] = read(joinpath(@__DIR__, "gh_probe.jl"), String)
    a["code_script"] = read(@__FILE__, String)
end

Random.seed!(20261109 + L)
t_start = now()
for (theta, g) in points
    grp = "L$(L)/theta$(theta)/g$(g)"
    h5open(fname, "cw") do f
        haskey(f, grp) && return
        SR = zeros(T + 1, NREAL); Isite = zeros(T + 1, L, NREAL); Iblock = zeros(T + 1, L, NREAL)
        ZZ = zeros(T + 1, NREAL); XX = similar(ZZ); YY = similar(ZZ)
        hs = zeros(L, NREAL); Js = zeros(L - 1, NREAL); cs = zeros(Int, L, NREAL)
        for n in 1:NREAL
            r = gh_probe_run(L, theta, g, T; x0 = x0)
            SR[:, n] = r.SR; Isite[:, :, n] = r.Isite; Iblock[:, :, n] = r.Iblock
            ZZ[:, n] = r.ZZ; XX[:, n] = r.XX; YY[:, n] = r.YY
            hs[:, n] = r.h; Js[:, n] = r.J; cs[:, n] = r.config
        end
        for (k, v) in pairs((; SR, Isite, Iblock, ZZ, XX, YY, h = hs, J = Js, config = cs))
            f["$grp/$k"] = v
        end
    end
    println("theta=$theta g=$g done, ", round(now() - t_start, Second))
end
