# Postprocessing of the canonical survey (script.canonical.jl): raw data on E: -> plotdata data/majorana/canonical_plot.hdf5.
# Per (L, point): (estimate, 2.5%, 97.5%) from 1000 bootstrap resamples of the realisations for
# eta = log10 median(P_e/P_o), gap = median scaled antipodal gap, edge/bulk = mean A_pi, chi = median chi_n.
# usage: julia post.canonical.jl [rawdir]   (default E:/info_lattice_fisher_information/rawdata/majorana)
using HDF5, Statistics, Random
Random.seed!(1)
raw = isempty(ARGS) ? "E:/info_lattice_fisher_information/rawdata/majorana" : ARGS[1]
out = joinpath(@__DIR__, "..", "..", "data", "majorana", "canonical_plot.hdf5")
function boot(stat, n; B = 1000)
    m = stat(1:n); s = [stat(rand(1:n, n)) for _ in 1:B]
    [m, quantile(s, 0.025), quantile(s, 0.975)]
end
h5open(out, "w") do g
    for L in (8, 10, 12)
        fn = joinpath(raw, "L$L", "canonical.hdf5"); isfile(fn) || continue
        h5open(fn, "r") do f
            for nm in keys(f["L$L"])
                grp = "L$L/$nm"
                e = read(f["$grp/psd_even"]); o = read(f["$grp/psd_odd"]); n = length(e); N = 2^L
                gap = reshape(read(f["$grp/gap"]), N, n); chi = reshape(read(f["$grp/chi"]), N, n)
                ae = read(f["$grp/edge_pi"]); ab = read(f["$grp/bulk_pi"])
                g["$grp/eta"] = boot(r -> log10(median(e[r] ./ o[r])), n)
                g["$grp/gap"] = boot(r -> median(gap[:, r]), n)
                g["$grp/edge"] = boot(r -> mean(ae[r]), n); g["$grp/bulk"] = boot(r -> mean(ab[r]), n)
                g["$grp/chi"] = boot(r -> median(chi[:, r]), n); g["$grp/n"] = n
            end
        end
    end
end
println("wrote ", out)
