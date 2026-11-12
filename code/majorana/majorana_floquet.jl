#####################################################################
# The DTC circuit at theta=0, h=0 (no gate noise) as a free Majorana
# Floquet chain (K.A.R.E.N., Ch5).
#
# Jordan-Wigner: g_{2j-1} = (prod_{k<j} X_k) Z_j,  g_{2j} = (prod_{k<j} X_k) Y_j
#   X_j        = i g_{2j-1} g_{2j}
#   Z_j Z_{j+1} = i g_{2j} g_{2j+1}
#   Z_1        = g_1   (edge spin = edge Majorana)
# A gate exp(phi g_a g_b) acts on the Majoranas (Heisenberg picture) as a
# rotation by 2 phi in the (a,b) plane:  g_a -> cos2phi g_a + sin2phi g_b,
# g_b -> cos2phi g_b - sin2phi g_a.
#   kick  exp(-i pi g/2 X_j)  = exp(+(pi g/2) g_{2j-1} g_{2j}) : rotation by  pi g
#   bond  exp(-i J_j Z_jZ_{j+1}) = exp(+J_j g_{2j} g_{2j+1})  : rotation by 2 J_j
# U_F = K D (D first), Heisenberg map U^dag g_a U = sum_b M_ab g_b, M = M_K M_D.
# Infinite temperature edge autocorrelator <Z_1(t) Z_1> = (M^t)_{11}.
#####################################################################

using LinearAlgebra

# Heisenberg matrix of exp(phi g_a g_b) on n Majoranas
function rot!(M::AbstractMatrix, a::Int, b::Int, phi::Real)
    c, s = cos(2phi), sin(2phi)
    R = Matrix{Float64}(I, size(M)...)
    R[a, a] = c; R[a, b] = s
    R[b, a] = -s; R[b, b] = c
    return R
end

"""
    floquet_majorana(g, J; pbc=0)

2L x 2L Heisenberg map M of U_F = exp(-i pi g/2 sum X) exp(-i sum J_j Z_jZ_{j+1}).
J has L-1 entries (open chain) or L (ring); for the ring the closing bond is
s * i g_{2L} g_1 with s = pbc = +1 or -1 (the fermionic boundary condition).
"""
function floquet_majorana(g::Real, J::AbstractVector; pbc::Int = 0)
    L = pbc == 0 ? length(J) + 1 : length(J)
    n = 2L
    MK = Matrix{Float64}(I, n, n)
    for j in 1:L
        a, b = 2j - 1, 2j
        c, s = cos(pi * g), sin(pi * g)
        MK[a, a] = c; MK[a, b] = s; MK[b, a] = -s; MK[b, b] = c
    end
    MD = Matrix{Float64}(I, n, n)
    for j in 1:L-1
        a, b = 2j, 2j + 1
        c, s = cos(2J[j]), sin(2J[j])
        MD[a, a] = c; MD[a, b] = s; MD[b, a] = -s; MD[b, b] = c
    end
    if pbc != 0
        # closing bond pbc * J_L * i g_{2L} g_1  -> exp(pbc J_L g_{2L} g_1)
        a, b = 2L, 1
        c, s = cos(2 * pbc * J[L]), sin(2 * pbc * J[L])
        MD[a, a] = c; MD[a, b] = s; MD[b, a] = -s; MD[b, b] = c
    end
    return MK * MD
end

# quasienergies in (-pi, pi] of the orthogonal map
quasienergies(M) = angle.(eigvals(M))

"""
    edge_modes(M; tol)

Number of eigenvalues of M at +1 (0 modes) and -1 (pi modes) within tol,
the edge weight (first and last 10% of Majoranas) of the pi mode eigenvectors.
"""
function edge_modes(M; tol = 1e-8)
    F = eigen(M)
    e = F.values
    n0 = count(x -> abs(x - 1) < tol, e)
    npi = count(x -> abs(x + 1) < tol, e)
    return (; n0, npi)
end

"""
    majorana_number(A)

Kitaev's Majorana number of a real antisymmetric single-particle generator on a
ring, from the two fermionic boundary conditions: sgn Pf(A_+) sgn Pf(A_-).
Here we use sign(Pf A) = sign of the product over 2x2 blocks of the real Schur form,
computed via the determinant-free route: Pf(A)^2 = det(A), sign from a skew LTL.
"""
function pfaffian_sign(A::AbstractMatrix)
    # Parlett-Reid style skew-symmetric elimination (real), returns sign of Pf
    A = copy(Matrix{Float64}(A))
    n = size(A, 1)
    sgn = 1.0
    for k in 1:2:n-1
        # pivot: largest |A[k, j]| for j > k
        jmax = k + argmax(abs.(A[k, k+1:n]))
        if jmax != k + 1
            # swap rows/cols k+1 and jmax
            A[[k + 1, jmax], :] = A[[jmax, k + 1], :]
            A[:, [k + 1, jmax]] = A[:, [jmax, k + 1]]
            sgn = -sgn
        end
        piv = A[k, k+1]
        piv == 0 && return 0.0
        sgn *= sign(piv)
        if k + 2 <= n
            tau = A[k, k+2:n] ./ piv
            # eliminate (Schur complement): A[k+2:n, k+2:n] -= tau * A[k+1, k+2:n]' - A[k+1, k+2:n] * tau'
            A[k+2:n, k+2:n] .-= tau * A[k+1, k+2:n]' .- A[k+1, k+2:n] * tau'
        end
    end
    return sgn
end

# real antisymmetric generator with branch cut at quasienergy pi (cut=:pi) or 0 (cut=:zero)
function generator(M; cut = :pi)
    X = cut == :pi ? M : -M
    A = real(log(Matrix{ComplexF64}(X)))
    return (A - A') / 2
end

"""
    floquet_invariants(g, J)

(nu_0, nu_pi) in {0,1}: (-1)^nu_0 = Majorana number of log(M_ring) (cut at pi),
(-1)^nu_pi = Majorana number of log(-M_ring) (cut at 0), J of length L (ring).
"""
function floquet_invariants(g, J)
    Mp = floquet_majorana(g, J; pbc = +1)
    Mm = floquet_majorana(g, J; pbc = -1)
    m0 = pfaffian_sign(generator(Mp; cut = :pi)) * pfaffian_sign(generator(Mm; cut = :pi))
    mpi = pfaffian_sign(generator(Mp; cut = :zero)) * pfaffian_sign(generator(Mm; cut = :zero))
    return (nu0 = m0 < 0 ? 1 : 0, nupi = mpi < 0 ? 1 : 0)
end

edge_autocorrelation(M, T) = (v = zeros(size(M, 1)); v[1] = 1.0; [((M')^t * v)[1] for t in 0:T])

#####################################################################
# Chiral (BDI) invariants at theta=0, h=0
# Gamma = +1 on odd Majoranas, -1 on even ones; every gate is a rotation
# between an odd and an even Majorana, so Gamma R Gamma = R^T. In the two
# symmetric time frames M1 = K^1/2 D K^1/2, M2 = D^1/2 K D^1/2 one has
# Gamma M Gamma = M^T, so A = log M anticommutes with Gamma and is off
# diagonal, A = [[0, h], [-h^T, 0]] in (odd, even). nu = winding of h,
# computed in real space (Mondragon-Shem, Song, Hughes):
#   nu = (1/|C|) sum_{j in C} (h^{-1} [X, h])_{jj}, C the central half of the ring.
# Then nu_0 = (nu_1 + nu_2)/2, nu_pi = (nu_1 - nu_2)/2 (Asboth and Obuse).
#####################################################################

function kick_map(g, L)
    MK = Matrix{Float64}(I, 2L, 2L)
    for j in 1:L
        a, b = 2j - 1, 2j; c, s = cos(pi * g), sin(pi * g)
        MK[a, a] = c; MK[a, b] = s; MK[b, a] = -s; MK[b, b] = c
    end
    MK
end
function bond_map(J, L; pbc = 1)
    MD = Matrix{Float64}(I, 2L, 2L)
    for j in 1:L
        j == L && pbc == 0 && break
        a, b = 2j, j == L ? 1 : 2j + 1
        Jj = j == L ? pbc * J[L] : J[j]
        c, s = cos(2Jj), sin(2Jj)
        MD[a, a] = c; MD[a, b] = s; MD[b, a] = -s; MD[b, b] = c
    end
    MD
end

function winding_realspace(A, L)
    odd = 1:2:2L; even = 2:2:2L
    h = A[odd, even]
    X = Diagonal(collect(1.0:L))
    W = h \ (X * h - h * X)
    C = (L ÷ 4 + 1):(3L ÷ 4)
    return sum(W[j, j] for j in C) / length(C)
end

function chiral_invariants(g, J)
    L = length(J)
    Kh = kick_map(g / 2, L); Dh = bond_map(J ./ 2, L)
    M1 = Kh * Dh * Kh          # Heisenberg maps compose in the same order as the gates here
    M2 = Dh * kick_map(g, L) * Dh
    A1 = real(log(Matrix{ComplexF64}(M1))); A1 = (A1 - A1') / 2
    A2 = real(log(Matrix{ComplexF64}(M2))); A2 = (A2 - A2') / 2
    n1 = winding_realspace(A1, L); n2 = winding_realspace(A2, L)
    return (nu1 = n1, nu2 = n2, nu0 = (n1 + n2) / 2, nupi = (n1 - n2) / 2)
end
