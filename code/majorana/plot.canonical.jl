# Fig. fig:majorana_converse: the canonical driven Kitaev chain in spin language (script.canonical.jl -> post.canonical.jl).
# Points on the x axis, colour by L (edge map).
# a. line split eta (eq:line_split), b. median scaled antipodal gap (eq:band_gap) with pi ln2,
# c. late time pi signal of the edge spin (filled) and of a bulk spin (open), d. eigenstate spin glass order chi.
using HDF5, CairoMakie
include(joinpath(@__DIR__, "canonical.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "mbl_dtc", "code", "fss_style.jl"))
fss_theme!()
f = h5open(joinpath(@__DIR__, "..", "..", "data", "majorana", "canonical_plot.hdf5"), "r")   # plotdata from post.canonical.jl
Ls = [L for L in (8, 10, 12) if haskey(f, "L$L")]
names = [p[1] for p in POINTS]
ticks = ["0 SG", "PM", "π SG", "π SG\nstrong", "0π PM", "π SG\nw≠Δ", "cluster\nπ", "cluster π\nstrong"]
cols = fss_colors(:edge, 3)
fig = Figure(size = (1100, 650))
axs = [Axis(fig[1, 1], ylabel = L"\eta"), Axis(fig[1, 2], ylabel = L"\tilde{\Delta}", yscale = log10),
       Axis(fig[2, 1], ylabel = L"A_\pi"), Axis(fig[2, 2], ylabel = L"\mathrm{median}\ \chi")]
for ax in axs
    ax.xticks = (1:length(names), ticks); ax.xticklabelsize = 11
end
summary = Dict()
for (k, L) in enumerate(Ls), (ip, nm) in enumerate(names)
    grp = "L$L/$nm"; haskey(f, grp) || continue
    st = (eta = read(f["$grp/eta"]), gap = read(f["$grp/gap"]), edge = read(f["$grp/edge"]), bulk = read(f["$grp/bulk"]), chi = read(f["$grp/chi"]))
    summary[(L, nm)] = st
    x = ip + 0.18 * (k - 2)
    for (ax, s, mk) in ((axs[1], st.eta, :circle), (axs[2], st.gap, :circle), (axs[3], st.edge, :circle), (axs[4], st.chi, :circle))
        errorbars!(ax, [x], [s[1]], [s[1] - s[2]], [s[3] - s[1]], color = cols[k])
        scatter!(ax, [x], [s[1]], color = cols[k], marker = mk, markersize = 10, label = ax === axs[1] ? L"L=%$L" : nothing)
    end
    errorbars!(axs[3], [x], [st.bulk[1]], [st.bulk[1] - st.bulk[2]], [st.bulk[3] - st.bulk[1]], color = cols[k])
    scatter!(axs[3], [x], [st.bulk[1]], color = :white, strokecolor = cols[k], strokewidth = 2, marker = :diamond, markersize = 10)
end
hlines!(axs[1], [0.0], color = :gray, linestyle = :dash)
hlines!(axs[2], [pi * log(2)], color = :gray, linestyle = :dash)
ylims!(axs[3], -0.05, 1.05); ylims!(axs[4], -0.02, 0.9)
Legend(fig[1:2, 3], [[MarkerElement(color = cols[k], marker = :circle, markersize = 10) for k in eachindex(Ls)];
                     MarkerElement(color = :gray, marker = :circle, markersize = 10);
                     MarkerElement(color = :white, strokecolor = :gray, strokewidth = 2, marker = :diamond, markersize = 10)],
       [[L"L=%$L" for L in Ls]; "edge spin"; "bulk spin"], framevisible = true, labelsize = 12)
for (i, pos) in enumerate(((1, 1), (1, 2), (2, 1), (2, 2)))
    panel_label!(fig[pos...], letter(i))
end
save(joinpath(@__DIR__, "..", "..", "..", "Thesis", "tex", "chapters", "Chapter_5", "plots", "majorana_converse.png"), fig, px_per_unit = 2)
for L in Ls, nm in names
    haskey(summary, (L, nm)) || continue
    s = summary[(L, nm)]; r(x) = round.(x, sigdigits = 3)
    println("L$L $(rpad(nm, 18)) n=$(read(f["L$L/$nm/n"])) eta $(r(s.eta)) gap $(r(s.gap)) edge $(r(s.edge)) bulk $(r(s.bulk)) chi $(r(s.chi))")
end
close(f)
