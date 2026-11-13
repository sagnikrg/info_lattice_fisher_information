# Gullans-Huse probe on the pi mode (script.gh_edge.jl): reference qubit R in a Bell pair with the edge qubit (filled)
# or a bulk qubit x0=L/2 (open), late times 100..200. Colour by L (edge map).
# a. pi signal of the joint R-qubit correlator, mean (-1)^t <Z_R Z_x0>; b. 0 signal mean <Z_R Z_x0>; c. I(R:x0) in bits.
using HDF5, Statistics, Random, CairoMakie
include(joinpath(@__DIR__, "..", "..", "..", "mbl_dtc", "code", "fss_style.jl"))
fss_theme!(); Random.seed!(1)
f = h5open(joinpath(@__DIR__, "..", "..", "data", "majorana", "gh_edge.hdf5"), "r")
Ls = [L for L in (8, 10) if haskey(f, "L$L")]
names = ["dtc", "dtc_theta", "pm", "scrambling", "SG0", "PM", "SGpi", "SGpi_strong", "PM0pi", "cluster_pi"]
ticks = ["DTC\nθ=0", "DTC\nθ=0.2", "PM", "scram-\nbling", "K: 0 SG", "K: PM", "K: π SG", "K: π SG\nstrong", "K: 0π PM", "K: cluster\nπ"]
cols = fss_colors(:edge, 3)
boot(x; B = 1000) = (m = mean(x); s = [mean(rand(x, length(x))) for _ in 1:B]; (m, quantile(s, 0.025), quantile(s, 0.975)))
fig = Figure(size = (1100, 820))
axs = [Axis(fig[i, 1], ylabel = y) for (i, y) in enumerate((L"\overline{(-1)^t\langle Z_RZ_{x_0}\rangle}", L"\overline{\langle Z_RZ_{x_0}\rangle}", L"\overline{I(R:x_0)}\ \mathrm{[bits]}"))]
for ax in axs; ax.xticks = (1:length(names), ticks); ax.xticklabelsize = 11; end
for (k, L) in enumerate(Ls), (ip, nm) in enumerate(names)
    haskey(f, "L$L/$nm") || continue
    for (side, dx, filled) in (("edge", -0.1, true), ("bulk", 0.1, false))
        x = ip + dx + 0.06 * (k - 1.5)
        for (ax, key) in zip(axs, ("Api_ZZ", "A0_ZZ", "I_R"))
            s = boot(read(f["L$L/$nm/$side/$key"]))
            errorbars!(ax, [x], [s[1]], [s[1] - s[2]], [s[3] - s[1]], color = cols[k])
            scatter!(ax, [x], [s[1]], color = filled ? cols[k] : :white, strokecolor = cols[k], strokewidth = 2,
                     marker = filled ? :circle : :diamond, markersize = 10)
        end
    end
end
for ax in axs[1:2]; hlines!(ax, [0.0], color = :gray, linestyle = :dash); ylims!(ax, -0.1, 1.08); end
ylims!(axs[3], -0.05, 2.1); vlines!.(axs, 4.5, color = :gray, linewidth = 1)
Legend(fig[1:3, 2], [[MarkerElement(color = cols[k], marker = :circle, markersize = 10) for k in eachindex(Ls)];
                     MarkerElement(color = :gray, marker = :circle, markersize = 10);
                     MarkerElement(color = :white, strokecolor = :gray, strokewidth = 2, marker = :diamond, markersize = 10)],
       [[L"L=%$L" for L in Ls]; "R on the edge qubit"; "R on a bulk qubit"], framevisible = true, labelsize = 12)
for i in 1:3; panel_label!(fig[i, 1], letter(i)); end
out = get(ENV, "GH_PLOT_OUT", joinpath(@__DIR__, "..", "..", "data", "majorana", "gh_edge_preview.png"))
save(out, fig, px_per_unit = 2)
for L in Ls, nm in names
    haskey(f, "L$L/$nm") || continue
    r(side, key) = round.(boot(read(f["L$L/$nm/$side/$key"])), digits = 3)
    println("L$L $(rpad(nm, 12)) edge piZZ $(r("edge","Api_ZZ")) 0ZZ $(r("edge","A0_ZZ")) I $(r("edge","I_R")) | bulk piZZ $(r("bulk","Api_ZZ")) 0ZZ $(r("bulk","A0_ZZ")) I $(r("bulk","I_R"))")
end
close(f)
