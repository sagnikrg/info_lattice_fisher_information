#####################################################################
# The thesis circuit (brickwall_dtc of the header) with the gate noise
# switchable, for exact diagonalisation checks of the Majorana picture.
# Same layer order as brickwall_dtc / brickwall: U = U_even U_odd,
# odd gates: (ZZ+theta) . fields, boundary kicks on sites 1 and L;
# even gates: kicks on both sites . (ZZ+theta).
# noise=true reproduces brickwall_dtc (theta + pi/50 N(0,1), Z noise pi/50 N(0,1)).
#####################################################################

using LinearAlgebra
include("../.header/.src/circuit/brickwall.jl")

function brickwall_dtc_ctrl(L::Int, theta::Float64, g::Float64, h::AbstractVector, J::AbstractVector; noise::Bool = false)
    fonez = kron(Z, Z); ftwox = kron(X, X); ftwoy = kron(Y, Y)
    th = noise ? theta + randn() * pi / 50 : theta
    FU = fill(fill(0.0 * im, 4, 4), L - 1)
    for j in 1:L-1
        delh = noise ? randn(4) * pi / 50 : zeros(4)
        int1 = kron(RZ(delh[1]), RZ(delh[2]))
        int2 = exp(-im * J[j] * fonez - im * th / 2 * (ftwox + ftwoy))
        int3 = kron(RZ(delh[3]), RZ(delh[4]))
        FU[j] = int3 * int2 * int1
    end
    ZRow = RZ.(h)
    for i in 1:2:L-1
        FU[i] = FU[i] * kron(ZRow[i], ZRow[i+1])
    end
    gk = pi * g
    for i in 2:2:L-1
        FU[i] = kron(RX(gk), RX(gk)) * FU[i]
    end
    FU[1] = kron(RX(gk), I(2)) * FU[1]
    FU[end] = kron(I(2), RX(gk)) * FU[end]
    return FU
end

circuit_ctrl(L, theta, g, h, J; noise = false) = brickwall(brickwall_dtc_ctrl(L, theta, g, h, J; noise = noise))

# Z on site k of L (kron order: site 1 most significant)
function site_op(P, k, L)
    return kron(Matrix{ComplexF64}(I, 2^(k - 1), 2^(k - 1)), P, Matrix{ComplexF64}(I, 2^(L - k), 2^(L - k)))
end

"""
    edge_autocorr_ed(U, k, T)

Infinite temperature <Z_k(t) Z_k>/2^L for t = 0..T from the dense Floquet unitary.
"""
function edge_autocorr_ed(U, k, T, L)
    Zk = site_op(ComplexF64[1 0; 0 -1], k, L)
    F = schur(Matrix{ComplexF64}(U))   # unitary Q even for degenerate eigenvalues (eigen() is not orthonormal there)
    V = F.Z; ph = diag(F.T)
    Zt = V' * Zk * V               # Z in the eigenbasis
    W = abs2.(Zt)                  # |<m|Z|n>|^2
    N = size(U, 1)
    return [real(sum(W .* (conj.(ph) .^ t * transpose(ph .^ t)))) / N for t in 0:T]
end
