#####################################################################
# "Are all Majorana pi modes time crystals?"
# The canonical driven (Floquet) Kitaev chain, translated back to spins
# with the Jordan-Wigner map of majorana_floquet.jl:
#   chemical potential  mu c^dag c      <->  X_j                (kick)
#   hopping + pairing  (w, Delta)       <->  (w+D)/2 ZZ + (w-D)/2 YY
#   next nearest       i g_{2j} g_{2j+3} <->  Z_j X_{j+1} Z_{j+2}  (cluster)
# U = K . B, K = prod_j exp(-i pi g_j/2 X_j), B = bonds (ZZ+YY) then cluster.
# Dense, open chain, site 1 most significant (as site_op in ed_circuit.jl).
# Diagnostics of Chapter 3 on the spins: PSD (even/odd lines), antipodal
# gap, plus edge/bulk Z autocorrelators and the eigenstate spin glass order.
#####################################################################

using LinearAlgebra, Statistics, Random

const sX = ComplexF64[0 1; 1 0]
const sY = ComplexF64[0 -im; im 0]
const sZ = ComplexF64[1 0; 0 -1]

# apply a k-site gate G on sites j..j+k-1 to all columns of Psi (in place copy)
function apply_gate(Psi::Matrix{ComplexF64}, G::Matrix{ComplexF64}, j::Int, k::Int, L::Int)
    N = size(Psi, 2)
    a = 2^(L - j - k + 1); b = 2^k; c = 2^(j - 1)
    T = reshape(Psi, a, b, c * N)
    out = similar(T)
    @inbounds for m in 1:c*N
        @views out[:, :, m] = T[:, :, m] * transpose(G)
    end
    return reshape(out, 2^L, N)
end

"""
    floquet_kitaev(L; g, Jz, Jy, Kc)

Spin form of the driven Kitaev chain. g::Vector (L), Jz, Jy::Vector (L-1), Kc::Vector (L-2).
"""
function floquet_kitaev(L; g, Jz, Jy = zeros(L - 1), Kc = zeros(L - 2))
    U = Matrix{ComplexF64}(I, 2^L, 2^L)
    for j in 1:L-1
        (Jz[j] == 0 && Jy[j] == 0) && continue
        G = exp(-im * (Jz[j] * kron(sZ, sZ) + Jy[j] * kron(sY, sY)))
        U = apply_gate(U, G, j, 2, L)
    end
    for j in 1:L-2
        Kc[j] == 0 && continue
        G = exp(-im * Kc[j] * kron(sZ, kron(sX, sZ)))
        U = apply_gate(U, G, j, 3, L)
    end
    for j in 1:L
        U = apply_gate(U, exp(-im * pi * g[j] / 2 * sX), j, 1, L)
    end
    return U
end

zdiag(k, L) = [iseven((s - 1) >> (L - k)) ? 1.0 : -1.0 for s in 1:2^L]   # Z_k on basis states

function psd_lines(om::Vector{Float64})
    N = length(om); w = sort(mod.(om, 2pi))
    f = N .* w ./ (2pi) .- (1:N)
    P = abs2.(fft_naive(f)) ./ N
    ks = (N ÷ 20):(N ÷ 2)   # window of the line split eta, eq:line_split
    ev = [P[k+1] for k in ks if iseven(k)]; od = [P[k+1] for k in ks if isodd(k)]
    return sum(ev), sum(od)
end
fft_naive(f) = [sum(f[n] * cis(-2pi * k * (n - 1) / length(f)) for n in eachindex(f)) for k in 0:length(f)-1]

function antipodal(om::Vector{Float64})
    N = length(om); w = sort(mod.(om, 2pi))
    d = similar(w)
    for n in 1:N
        x = mod(w[n] + pi, 2pi)
        i = searchsortedfirst(w, x)
        c1 = w[mod1(i, N)]; c2 = w[mod1(i - 1, N)]
        d[n] = min(abs(mod(c1 - x + pi, 2pi) - pi), abs(mod(c2 - x + pi, 2pi) - pi))
    end
    return N .* d
end

# window average over t = t1..t2 of exp(i w t)
kern(w, t1, t2) = abs(mod(w + pi, 2pi) - pi) < 1e-12 ? 1.0 + 0im :
                  (cis(w * t1) - cis(w * (t2 + 1))) / ((1 - cis(w)) * (t2 - t1 + 1))

function diagnose(U, L; t1 = 200, t2 = 400)
    F = schur(U); V = F.Z; om = angle.(diag(F.T)); N = 2^L
    ev, od = psd_lines(om)
    gap = antipodal(om)
    dw = [om[n] - om[m] for m in 1:N, n in 1:N]                  # e^{i(om_n - om_m) t}
    Kpi = real.(kern.(dw .+ pi, t1, t2)); K0 = real.(kern.(dw, t1, t2))
    A = Dict{Symbol,Float64}()
    for (nm, k) in ((:edge, 1), (:bulk, L ÷ 2))
        z = zdiag(k, L)
        W = abs2.(V' * (z .* V))
        A[Symbol(nm, :_pi)] = sum(W .* Kpi) / N
        A[Symbol(nm, :_0)] = sum(W .* K0) / N
    end
    # eigenstate spin glass order chi = sum_{i != j} <Z_i Z_j>^2 / (L(L-1))
    Pr = abs2.(V); Zs = hcat([zdiag(k, L) for k in 1:L]...)
    chi = zeros(N)
    for i in 1:L, j in i+1:L
        chi .+= 2 .* (transpose(Zs[:, i] .* Zs[:, j]) * Pr)[:] .^ 2
    end
    chi ./= L * (L - 1)
    return (psd_even = ev, psd_odd = od, gap = gap, chi = chi, A...)
end

# the canonical points (name, label, g, Jz, Jy, Kc, relative disorder w of the couplings);
# couplings x(1 +- w), kick g +- 0.025
const POINTS = [
    (:SG0,        "0 mode, spin glass",              0.1, pi / 4,  0.0,     0.0,    0.25),
    (:PM,         "trivial paramagnet",              0.5, pi / 8,  0.0,     0.0,    0.25),
    (:SGpi,       "pi mode, spin glass (DTC)",       0.9, pi / 4,  0.0,     0.0,    0.25),
    (:SGpi_strong,"pi mode, spin glass, strong dis.",0.9, pi / 4,  0.0,     0.0,    0.8),
    (:PM0pi,      "0 and pi modes, paramagnet",      0.5, 3pi / 8, 0.0,     0.0,    0.25),
    (:SGpi_aniso, "pi mode, w != Delta",             0.9, pi / 4,  pi / 32, 0.0,    0.25),
    (:cluster_pi, "two pi modes, cluster",           0.9, 0.0,     0.0,     pi / 4, 0.25),
    (:cluster_pi_strong, "two pi modes, cluster, strong dis.", 0.9, 0.0, 0.0, pi / 4, 0.8),
]

function realisation(pt, L, rng)
    _, _, g, Jz, Jy, Kc, w = pt
    gv = g .+ 0.05 .* (rand(rng, L) .- 0.5)
    d(x, n) = x .* (1 .+ 2w .* (rand(rng, n) .- 0.5))
    return floquet_kitaev(L; g = gv, Jz = d(Jz, L - 1), Jy = d(Jy, L - 1), Kc = d(Kc, L - 2))
end
