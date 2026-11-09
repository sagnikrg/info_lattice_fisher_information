#####################################################################
# Gullans-Huse type reference qubit probe for the DTC Floquet circuit
#
# Gullans and Huse, PRL 125, 070606 (2020): a reference qubit R starts
# maximally entangled with one site x0 of the chain, and the dynamics
# acts on the chain only. Here the dynamics is our Floquet circuit
# (circuit_dtc from the header, the code defines the model), without
# measurements. Since U acts on the chain alone, S(R)=1 bit at all times,
# so the information sits in how R is correlated with parts of the chain:
# the mutual informations I(R:A) and the R-site correlators below.
#
# Layout of the state: a 2 x 2^L matrix M[r, s], r the reference, s the
# chain in the kron order of circuit_dtc (site 1 most significant).
# One Floquet period: M -> M * transpose(U).
#####################################################################

using LinearAlgebra

include("../.header/.src/circuit/brickwall.jl")
include("../.header/.src/circuit/time_crystals.jl")

# bit of site k (1..L) in the 0-based chain index s
site_bit(s, k, L) = (s >> (L - k)) & 1

"""
    gh_initial_state(L, x0, config)

Reference R in the Bell pair (|0>_R|0>_x0 + |1>_R|1>_x0)/sqrt(2) with site x0,
every other site j in the Z eigenstate config[j] (0 or 1).
"""
function gh_initial_state(L::Int, x0::Int, config::AbstractVector{<:Integer})
    M = zeros(ComplexF64, 2, 2^L)
    for r in 0:1
        s = 0
        for k in 1:L
            b = k == x0 ? r : config[k]
            s |= b << (L - k)
        end
        M[r+1, s+1] = 1 / sqrt(2)
    end
    return M
end

# tensor axis of site k when vec(M) is reshaped to (2,2,...,2): axis 1 is R,
# and the chain index is column major, so site L is axis 2 and site 1 is axis L+1
site_axis(k, L) = L + 2 - k

"""
    reduced_dm(M, keep_R, sites, L)

Reduced density matrix of R (if keep_R) and the chain sites `sites`.
"""
function reduced_dm(M::AbstractMatrix, keep_R::Bool, sites::AbstractVector{<:Integer}, L::Int)
    T = reshape(M, ntuple(_ -> 2, L + 1))
    keep = Int[]
    keep_R && push!(keep, 1)
    append!(keep, site_axis.(sites, L))
    rest = setdiff(1:L+1, keep)
    A = reshape(permutedims(T, (keep..., rest...)), 2^length(keep), :)
    return A * A'
end

function entropy_bits(rho::AbstractMatrix)
    p = eigvals(Hermitian(rho))
    return -sum(x -> x > 1e-14 ? x * log2(x) : 0.0, p)
end

"""
    entropy_part(M, keep_R, sites, L)

Entropy in bits of R (if keep_R) and `sites`, from the singular values of the
bipartition of the pure global state, i.e. always on the smaller side.
"""
function entropy_part(M::AbstractMatrix, keep_R::Bool, sites::AbstractVector{<:Integer}, L::Int)
    T = reshape(M, ntuple(_ -> 2, L + 1))
    keep = Int[]
    keep_R && push!(keep, 1)
    append!(keep, site_axis.(sites, L))
    (isempty(keep) || length(keep) == L + 1) && return 0.0
    rest = setdiff(1:L+1, keep)
    A = reshape(permutedims(T, (keep..., rest...)), 2^length(keep), :)
    p = svdvals(A) .^ 2
    return -sum(x -> x > 1e-14 ? x * log2(x) : 0.0, p)
end

"""
    mutual_info_R(M, sites, L)

I(R:A) = S(R) + S(A) - S(RA) in bits, A = `sites`. At most 2 bits.
"""
function mutual_info_R(M, sites, L)
    return entropy_part(M, true, Int[], L) + entropy_part(M, false, sites, L) -
           entropy_part(M, true, sites, L)
end

const PX = ComplexF64[0 1; 1 0]
const PY = ComplexF64[0 -im; im 0]
const PZ = ComplexF64[1 0; 0 -1]

"""
    R_site_correlator(M, P, k, L)

<P_R P_k> from the two site reduced density matrix of R and site k.
The kron order of reduced_dm with keep=(R, k) is k (x) R, R the fastest index.
"""
function R_site_correlator(M, P, k, L)
    rho = reduced_dm(M, true, [k], L)
    return real(tr(rho * kron(P, P)))
end

# block of w sites centred at x0, clipped to the chain
function centred_block(x0, w, L)
    a = clamp(x0 - (w - 1) ÷ 2, 1, L - w + 1)
    return collect(a:a+w-1)
end

"""
    gh_probe_run(L, theta, g, T; x0=L÷2, h, J, config)

One realisation: build U = circuit_dtc(L, theta, 1-g, h, J), evolve T periods and
record at every period t=0..T
  SR[t]          entropy of R (stays 1, a check)
  Isite[t, j]    I(R:j) for every site j
  Iblock[t, w]   I(R:block of w sites centred at x0), w = 1..L
  ZZ, XX, YY[t]  <P_R P_x0>
"""
function gh_probe_run(L::Int, theta::Float64, g::Float64, T::Int;
                      x0::Int = L ÷ 2,
                      h = rand(L) * 2pi, J = rand(L - 1) * pi,
                      config = rand(0:1, L))
    U = circuit_dtc(L, theta, 1 - g, h, J)
    Ut = transpose(U)
    M = gh_initial_state(L, x0, config)

    SR = zeros(T + 1)
    Isite = zeros(T + 1, L)
    Iblock = zeros(T + 1, L)
    ZZ = zeros(T + 1); XX = zeros(T + 1); YY = zeros(T + 1)

    for t in 0:T
        SR[t+1] = entropy_part(M, true, Int[], L)
        for j in 1:L
            Isite[t+1, j] = mutual_info_R(M, [j], L)
        end
        for w in 1:L
            Iblock[t+1, w] = mutual_info_R(M, centred_block(x0, w, L), L)
        end
        ZZ[t+1] = R_site_correlator(M, PZ, x0, L)
        XX[t+1] = R_site_correlator(M, PX, x0, L)
        YY[t+1] = R_site_correlator(M, PY, x0, L)
        t < T && (M = M * Ut)
    end
    return (; SR, Isite, Iblock, ZZ, XX, YY, h, J, config)
end
