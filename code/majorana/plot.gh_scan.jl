# Fig. fig:gh_probe_scan: the reference qubit probe along g for theta = 0..0.4 (script.gh_scan.jl).
# a. R on the edge qubit, b. R on the bulk qubit L/2: A_0^R for g < 1/2 and A_pi^R for g > 1/2 (late time window 100..200);
# c. I(R:x0) of the bulk qubit. Colours: theta (edge map), line style: L (dashed 8, solid 10). 95% bootstrap bands.
using HDF5, Statistics, Random, CairoMakie
include(joinpath(@__DIR__, "..", "..", "..", "mbl_dtc", "code", "fss_style.jl"))
fss_theme!(); Random.seed!(1)
f = h5open(joinpath(@__DIR__, "..", "..", "data", "majorana", "gh_scan.hdf5"), "r")
gs = read(f["g"]); Ls = [L for L in (8, 10, 12) if haskey(f, "L$L")]
ths = (0.0, 0.1, 0.2, 0.3, 0.4); cols = fss_colors(:edge, 6)[1:5]; styles = Dict(8 => :dash, 10 => :solid, 12 => :dot)
boot(x; B = 1000) = (m = mean(x); s = [mean(rand(x, length(x))) for _ in 1:B]; (m, quantile(s, 0.025), quantile(s, 0.975)))
fig = Figure(size = (1350, 400))
ylab = L"A^R_0\ (g<\frac{1}{2}),\quad A^R_\pi\ (g>\frac{1}{2})"
axs = [Axis(fig[1, 1], xlabel = L"g", ylabel = ylab, title = "R on the edge qubit"),
       Axis(fig[1, 2], xlabel = L"g", ylabel = ylab, title = "R on the bulk qubit L/2"),
       Axis(fig[1, 3], xlabel = L"g", ylabel = L"\overline{I(R:x_0)}\ \mathrm{[bits]}", title = "R on the bulk qubit L/2")]
for L in Ls, (it, th) in enumerate(ths)
    grp = "L$L/theta$th"; haskey(f, grp) || continue
    for (ia, side) in enumerate(("edge", "bulk"))
        A = ifelse.(gs .< 0.5, read(f["$grp/$side/A0_ZZ"]), read(f["$grp/$side/Api_ZZ"]))
        st = [boot(A[i, :]) for i in eachindex(gs)]
        for sel in (gs .< 0.5, gs .> 0.5)
            fss_lines!(axs[ia], gs[sel], [first.(st)[sel]], [getindex.(st, 2)[sel]], [getindex.(st, 3)[sel]], [cols[it]], [nothing];
                       linestyle = styles[L], markersize = 5, linewidth = L == 8 ? 2 : 3)
        end
    end
    I = read(f["$grp/bulk/I_R"]); st = [boot(I[i, :]) for i in eachindex(gs)]
    fss_lines!(axs[3], gs, [first.(st)], [getindex.(st, 2)], [getindex.(st, 3)], [cols[it]], [nothing];
               linestyle = styles[L], markersize = 5, linewidth = L == 8 ? 2 : 3)
end
for ax in axs[1:2]; ylims!(ax, -0.05, 1.08); vlines!(ax, [0.5], color = :gray, linestyle = :dash); end
ylims!(axs[3], -0.05, 2.05)
for i in 1:3; panel_label!(fig[1, i], letter(i)); end
Legend(fig[1, 4], [[LineElement(color = cols[i], linewidth = 3) for i in eachindex(ths)]; [LineElement(color = :gray, linestyle = styles[L], linewidth = 3) for L in Ls]],
       [[L"\theta=%$th" for th in ths]; [L"L=%$L" for L in Ls]], framevisible = true, labelsize = 12)
save(get(ENV, "GH_SCAN_OUT", joinpath(@__DIR__, "..", "..", "..", "Thesis", "tex", "chapters", "Chapter_5", "plots", "gh_probe_scan.png")), fig, px_per_unit = 2)
for L in Ls, th in ths
    grp = "L$L/theta$th"; haskey(f, grp) || continue
    r(k, s) = round.(vec(mean(read(f["$grp/$s/$k"]), dims = 2)), digits = 2)
    println(grp, "\n  edge pi ", r("Api_ZZ", "edge"), "\n  edge 0  ", r("A0_ZZ", "edge"), "\n  bulk pi ", r("Api_ZZ", "bulk"), "\n  bulk 0  ", r("A0_ZZ", "bulk"), "\n  bulk I  ", r("I_R", "bulk"))
end
close(f)
