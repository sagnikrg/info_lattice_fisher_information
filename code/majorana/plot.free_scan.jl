# Fig. fig:majorana_invariants: a. the winding numbers nu_0, |nu_pi| against g (median over realisations,
# band = 95% bootstrap of the median), b. the typical splitting of the 0 and pi edge pairs, log10, against g.
using HDF5, Statistics, Random, CairoMakie
include(joinpath(@__DIR__, "..", "..", "..", "mbl_dtc", "code", "fss_style.jl"))
fss_theme!()
f = h5open(joinpath(@__DIR__, "..", "..", "data", "majorana", "free_scan.hdf5"), "r")
gs = read(f["g"]); Ls = sort([parse(Int, k[2:end]) for k in keys(f) if startswith(k, "L")])
Random.seed!(1)
function boot(x, stat; B = 1000)
    x = filter(isfinite, x); isempty(x) && return (NaN, NaN, NaN)
    s = [stat(rand(x, length(x))) for _ in 1:B]
    return (stat(x), quantile(s, 0.025), quantile(s, 0.975))
end
cols = fss_colors(:nu, length(Ls)); cols2 = fss_colors(:edge, length(Ls))
fig = Figure(size = (1150, 380))
ax1 = Axis(fig[1, 1], xlabel = L"g", ylabel = L"\nu", xgridvisible = true, ygridvisible = true)
ax2 = Axis(fig[1, 2], xlabel = L"g", ylabel = L"\log_{10}\,\delta\varepsilon")
for (k, L) in enumerate(Ls)
    for (q, ls, tr) in (("nu0", :solid, identity), ("nupi", :dash, abs))
        d = tr.(read(f["L$L/$q"]))
        st = [boot(d[i, :], median) for i in eachindex(gs)]
        fss_lines!(ax1, gs, [first.(st)], [getindex.(st, 2)], [getindex.(st, 3)], [cols[k]],
                   [q == "nu0" ? L"\nu_0,\ L=%$L" : L"|\nu_\pi|,\ L=%$L"]; linestyle = ls, markersize = 5)
    end
    for (q, ls) in (("split0", :solid), ("splitpi", :dash))
        d = log10.(max.(read(f["L$L/$q"]), 1e-16))
        st = [boot(d[i, :], median) for i in eachindex(gs)]
        fss_lines!(ax2, gs, [first.(st)], [getindex.(st, 2)], [getindex.(st, 3)], [cols2[k]],
                   [q == "split0" ? L"0,\ L=%$L" : L"\pi,\ L=%$L"]; linestyle = ls, markersize = 5)
    end
end
ylims!(ax1, -0.1, 1.2)
Legend(fig[1, 3], ax1, framevisible = true, labelsize = 11)
Legend(fig[1, 4], ax2, framevisible = true, labelsize = 11)
panel_label!(fig[1, 1], "a."); panel_label!(fig[1, 2], "b.")
out = joinpath(@__DIR__, "..", "..", "..", "Thesis", "tex", "chapters", "Chapter_5", "plots", "majorana_invariants.png")
save(out, fig, px_per_unit = 2)
for (k, L) in enumerate(Ls)
    for q in ("nu0", "nupi")
        d = read(f["L$L/$q"]); println("L=$L $q median: ", round.([median(filter(isfinite, abs.(d[i, :]))) for i in eachindex(gs)], digits = 2))
    end
    println("L=$L log10 split pi median: ", round.([median(log10.(max.(read(f["L$L/splitpi"])[i, :], 1e-16))) for i in eachindex(gs)], digits = 1))
    println("L=$L log10 split 0 median: ", round.([median(log10.(max.(read(f["L$L/split0"])[i, :], 1e-16))) for i in eachindex(gs)], digits = 1))
end
close(f)
