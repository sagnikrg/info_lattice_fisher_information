# ED survey of the canonical driven Kitaev chain in spin language (canonical.jl)
# usage: julia -t 1 script.canonical.jl  -> raw E:/info_lattice_fisher_information/rawdata/majorana/L$L/canonical.hdf5,
#        log data/majorana/log_canonical.txt; then post.canonical.jl -> plotdata
using HDF5
include("canonical.jl")
const OUT = joinpath(@__DIR__, "../../data/majorana")
const RAW = "E:/info_lattice_fisher_information/rawdata/majorana"
const NR = Dict(8 => 200, 10 => 30, 12 => 5)
logf(s) = (open(joinpath(OUT, "log_canonical.txt"), "a") do io; println(io, s); end; println(s))
for L in (8, 10, 12), pt in POINTS
    rng = MersenneTwister(1000L + findfirst(==(pt), POINTS))
    t = @elapsed rs = [diagnose(realisation(pt, L, rng), L) for _ in 1:NR[L]]
    mkpath(joinpath(RAW, "L$L"))
    h5open(joinpath(RAW, "L$L", "canonical.hdf5"), "cw") do f
        gname = "L$L/$(pt[1])"; haskey(f, gname) && delete_object(f, gname)
        gr = create_group(f, gname)
        for k in (:psd_even, :psd_odd, :edge_pi, :edge_0, :bulk_pi, :bulk_0)
            gr[string(k)] = [r[k] for r in rs]
        end
        gr["gap"] = vcat([r.gap for r in rs]...); gr["chi"] = vcat([r.chi for r in rs]...)
    end
    pr = log10(median([r.psd_even / r.psd_odd for r in rs]))
    f(s) = round(mean(r[s] for r in rs), digits = 3)
    logf("L$L $(rpad(pt[1], 18)) eta $(round(pr, digits = 2)) gap $(round(median(vcat([r.gap for r in rs]...)), sigdigits = 3)) chi $(round(median(vcat([r.chi for r in rs]...)), digits = 3)) edge pi/0 $(f(:edge_pi))/$(f(:edge_0)) bulk pi/0 $(f(:bulk_pi))/$(f(:bulk_0))  $(round(t, digits = 1))s")
end
