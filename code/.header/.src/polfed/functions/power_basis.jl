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

function apply_polynomial(a::AbstractVector, k::Int, phi::Real)
    n=length(a)
    out=zeros(ComplexF64, n+k)
    for m in 0:k
        out[m+1:m+n] .+= exp(-im*m*phi) .* a
    end
    return out
end


#############
# Arnoldi in the power basis
###########

# returns the Hessenberg matrix H (ncv x ncv) and the Krylov vectors as columns of V (coefficients,
# zero padded to n = ncv*k+1), from the overlaps c (length >= n). S = overlap(c[1:n]) is kept as a
# Toeplitz matrix (S*v is an FFT product, never a dense n x n matrix), and S*v_i is kept for every
# Krylov vector, so <v_i|w> = (S v_i)' w costs a dot product only.

function arnoldi_power(c::AbstractVector, k::Int, ncv::Int, phi::Real)

    n=ncv*k+1
    V=zeros(ComplexF64, n, ncv)
    SV=zeros(ComplexF64, n, ncv)
    H=zeros(ComplexF64, ncv+1, ncv)
    Sn=overlap(c[1:n])

    V[1, 1]=1.0                                  # |v_1> = |E>
    SV[:, 1]=Sn*V[:, 1]

    for j in 1:ncv
        len=(j-1)*k+1                            # support of v_j
        w=zeros(ComplexF64, n)
        w[1:len+k]=apply_polynomial(V[1:len, j], k, phi)

        for pass in 1:2                          # classical Gram-Schmidt, twice
            h=SV[:, 1:j]'*w
            H[1:j, j].+=h
            w.-=V[:, 1:j]*h
        end

        H[j+1, j]=sqrt(max(real(dot(w, Sn*w)), 0.0))
        if j<ncv
            abs(H[j+1, j]) < 1e-14 && return H[1:j, 1:j], V[:, 1:j]
            V[:, j+1]=w ./ H[j+1, j]
            SV[:, j+1]=Sn*V[:, j+1]
        end
    end

    return H[1:ncv, 1:ncv], V
end
