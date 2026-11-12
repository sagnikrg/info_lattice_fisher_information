# Fig. fig:majorana_edge: late time edge signals from ED (noise off), L >= 8.
# Columns: h = 0 and h random; on the left of each panel the 0 mode signal A_0 (g < 1/2), on the right the pi mode signal A_pi.
# Colours: L (edge map), line style: theta.
using HDF5, Statistics, Random, CairoMakie
include(joinpath(@__DIR__, "..", "..", "..", "mbl_dtc", "code", "fss_style.jl"))
fss_theme!(); Random.seed!(1)
f = h5open(joinpath(@__DIR__, "..", "..", "data", "majorana", "edge_ed.hdf5"), "r")
gs = read(f["g"]); Ls = [L for L in (8, 10, 12) if haskey(f, "L$L")]
boot(x; B = 1000) = (m = mean(x); s = [mean(rand(x, length(x))) for _ in 1:B]; (m, quantile(s, 0.025), quantile(s, 0.975)))
cols = fss_colors(:edge, length(Ls)); styles = Dict(0.0 => :solid, 0.1 => :dash, 0.2 => :dot, 0.3 => :dashdot)
fig = Figure(size = (1000, 380))
for (ip, hon) in enumerate((0, 1))
    ax = Axis(fig[1, ip], xlabel = L"g", ylabel = L"A_0\ \ (g<\frac{1}{2}),\quad A_\pi\ \ (g>\frac{1}{2})")
    for (k, L) in enumerate(Ls), th in (0.0, 0.1, 0.2, 0.3)
        grp = "L$L/h$hon/theta$th"; haskey(f, grp) || continue
        A = ifelse.(gs .< 0.5, read(f["$grp/A0"]), read(f["$grp/Api"]))
        st = [boot(A[i, :]) for i in eachindex(gs)]
        for side in (gs .< 0.5, gs .> 0.5)
            fss_lines!(ax, gs[side], [first.(st)[side]], [getindex.(st, 2)[side]], [getindex.(st, 3)[side]], [cols[k]],
                       [side[1] ? L"L=%$L,\ \theta=%$th" : nothing]; linestyle = styles[th], markersize = 5)
        end
    end
    ylims!(ax, -0.05, 1.0); vlines!(ax, [0.5], color = :gray, linestyle = :dash)
    ip == 2 && Legend(fig[1, 3], ax, framevisible = true, labelsize = 11, merge = true, unique = true)
    panel_label!(fig[1, ip], letter(ip))
end
save(joinpath(@__DIR__, "..", "..", "..", "Thesis", "tex", "chapters", "Chapter_5", "plots", "majorana_edge.png"), fig, px_per_unit = 2)
for L in Ls, hon in (0, 1), th in (0.0, 0.1, 0.2, 0.3)
    grp = "L$L/h$hon/theta$th"; haskey(f, grp) || continue
    println(grp, " n=", size(read(f["$grp/A0"]), 2), " A0=", round.(vec(mean(read(f["$grp/A0"]), dims = 2)), digits = 3), " Api=", round.(vec(mean(read(f["$grp/Api"]), dims = 2)), digits = 3))
end
close(f)
