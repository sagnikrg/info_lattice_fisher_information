#######################################
# The power basis |j> = U^j|E>: overlaps, inner products, g_k(U) and Arnoldi on coefficient vectors
# (used by polfed_mps.jl; needs functions/overlap.jl, states.jl, operator.jl)
#######################################

# the idea of the power basis and the notation (k, ncv) are explained at the top of polfed_mps.jl


#####################################
# The first column of the overlap matrix: one time evolution
#####################################

function overlaps(U, psi0, N::Int)

    c=zeros(ComplexF64, N+1)
    psi=psi0
    ref=psi0                                     # psi0 as the same type as psi (dense once switched)

    c[1]=inner_state(ref, psi)
    for j in 1:N
        psi=floquet_step(U, psi)
        (psi isa ITensor && ref isa MPS) && (ref=prod(psi0))     # DenseSwitch: compare against the dense psi0
        c[j+1]=inner_state(ref, psi)
    end

    c./=real(c[1])
    c[1]=1.0                                     # c_0 = <E|E> = 1 exactly (Toeplitz needs it real)
    return c
end


#####################################
# Inner products in the power basis
#####################################

# <a|b> = sum_{i,j} a_i^* b_j <i|j>, coefficient vectors of any length (zero padded)
# overlap: the overlap matrix (at least as large as the vectors), sliced dense, so for small checks

function inner_prod(vi::AbstractVector, vj::AbstractVector, overlap)
    n=max(length(vi), length(vj))
    a=vcat(vi, zeros(ComplexF64, n-length(vi)))
    b=vcat(vj, zeros(ComplexF64, n-length(vj)))
    return dot(a, overlap[1:n, 1:n]*b)
end

# the same for rows i, j of a coefficient matrix (one Krylov vector per row)

inner_prod(i::Int, j::Int, coeff::AbstractMatrix, overlap)=inner_prod(coeff[i, :], coeff[j, :], overlap)


#############
# Applying the polynomial on a vector in the power basis
###########

# g_k(U) sum_i a_i U^i|E> = sum_i sum_{m=0}^{k} e^{-i m phi} a_i U^{i+m}|E>, i.e a shift and sum:
#   (g_k a)_n = sum_{m=0}^{k} e^{-i m phi} a_{n-m},   length(g_k a) = length(a) + k
# a sliding window, computed as a running sum in O(n) instead of k+1 shifted copies in O(k n):
#   s_n = e^{-i phi} s_{n-1} + a_n - e^{-i (k+1) phi} a_{n-k-1}

function apply_polynomial(a::AbstractVector, k::Int, phi::Real)
    n=length(a)
    out=zeros(ComplexF64, n+k)
    z=exp(-im*phi); zk=exp(-im*(k+1)*phi)
    s=zero(ComplexF64)
    for i in 1:n+k
        s=z*s
        i <= n && (s+=a[i])
        i > k+1 && (s-=zk*a[i-k-1])
        out[i]=s
    end
    return out
end


#############
# Arnoldi in the power basis
###########

# returns the Hessenberg matrix H (ncv x ncv) and the Krylov vectors as columns of V (coefficients,
# zero padded to n = ncv*k+1), from the overlaps c (length >= n). S = overlap(c[1:n]) is a Toeplitz
# matrix, factorised once (S*v as an FFT product with a cached plan, never a dense n x n matrix), and
# S*v_i is kept for every Krylov vector, so <v_i|w> = (S v_i)' w costs a dot product only.
# The j-th Krylov vector has (j-1)k+1 nonzero coefficients (staircase), so the Gram-Schmidt only runs
# over the support of the new vector, and S w is computed once per step (for its norm and for S v_{j+1}).

function arnoldi_power(c::AbstractVector, k::Int, ncv::Int, phi::Real)

    n=ncv*k+1
    V=zeros(ComplexF64, n, ncv)
    SV=zeros(ComplexF64, n, ncv)
    H=zeros(ComplexF64, ncv+1, ncv)
    F=overlap_factorized(c, n)

    V[1, 1]=1.0                                  # |v_1> = |E>
    toeplitz_mul!(view(SV, :, 1), F, view(V, :, 1))
    Sw=zeros(ComplexF64, n)

    for j in 1:ncv
        len=(j-1)*k+1                            # support of v_j
        lw=min(len+k, n)                         # support of g_k v_j, and of everything subtracted from it
        w=zeros(ComplexF64, n)
        w[1:lw]=apply_polynomial(view(V, 1:len, j), k, phi)[1:lw]

        for pass in 1:2                          # classical Gram-Schmidt, twice
            h=view(SV, 1:lw, 1:j)'*view(w, 1:lw)
            H[1:j, j].+=h
            view(w, 1:lw) .-= view(V, 1:lw, 1:j)*h
        end

        toeplitz_mul!(Sw, F, w)
        H[j+1, j]=sqrt(max(real(dot(view(w, 1:lw), view(Sw, 1:lw))), 0.0))
        if j<ncv
            abs(H[j+1, j]) < 1e-14 && return H[1:j, 1:j], V[:, 1:j]
            V[:, j+1]=w ./ H[j+1, j]
            SV[:, j+1]=Sw ./ H[j+1, j]
        end
    end

    return H[1:ncv, 1:ncv], V
end
