# Quick text summary of gh_probe_small_L<L>.hdf5: late time (t in 50..100) means over
# realisations of I(R:x0), I(R:block of w sites around x0), and the R-x0 correlators
# at even and odd periods (period doubling shows up as opposite signs).
using HDF5, Statistics, Printf
L = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 8
f = h5open(joinpath(@__DIR__, "..", "..", "data", "gh_probe", "gh_probe_small_L$(L).hdf5"), "r")
x0 = read(HDF5.attributes(f)["x0"])
g0 = f["L$(L)"]
rows = []
for th in keys(g0), gg in keys(g0[th])
    d = g0[th][gg]
    theta = parse(Float64, th[6:end]); g = parse(Float64, gg[2:end])
    Is = read(d["Isite"]); Ib = read(d["Iblock"]); ZZ = read(d["ZZ"]); XX = read(d["XX"])
    late = 51:101
    ev = late[iseven.(late .- 1)]; od = late[isodd.(late .- 1)]
    push!(rows, (theta, g, mean(Is[late, x0, :]), [mean(Ib[late, w, :]) for w in 1:L],
                 mean(ZZ[ev, :]), mean(ZZ[od, :]), mean(XX[ev, :]), mean(XX[od, :]), size(Is, 3)))
end
sort!(rows, by = r -> (r[1], r[2]))
println("L=$L, x0=$x0, late times t=50..100, means over realisations")
@printf("%6s %5s %4s | %7s | %s | %7s %7s | %7s %7s\n", "theta", "g", "n", "I(R:x0)", join([@sprintf("Ib w=%-2d", w) for w in 1:L], " "), "ZZ even", "ZZ odd", "XX even", "XX odd")
for r in rows
    @printf("%6.2f %5.2f %4d | %7.3f | %s | %7.3f %7.3f | %7.3f %7.3f\n", r[1], r[2], r[9], r[3], join([@sprintf("%7.3f", x) for x in r[4]], " "), r[5], r[6], r[7], r[8])
end
close(f)
