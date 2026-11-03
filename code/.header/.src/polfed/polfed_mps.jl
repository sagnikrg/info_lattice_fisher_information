#######################################
# POLFED in the power basis (POLFED with matrix product states, Chapter 7)
#######################################

###
# Every Krylov vector of the Arnoldi algorithm on g_k(U) is a linear combination of
#
#   |0> = |E>,  |1> = U|E>,  |2> = U^2|E>, ...,  |N> = U^N|E>
#
# where |E> is the (random) initial state. Working in this (non-orthogonal, overcomplete) basis,
# a vector is just its coefficients, and everything the Arnoldi algorithm needs reduces to
#   - inner products, through the overlap matrix <i|j> = c_{j-i}, c_j = <E|U^j|E> (Hermitian Toeplitz)
#   - applying g_k(U), a shift and sum of the coefficients
#   - U|v>, a shift by one, for the eigenvalues <v|U|v>/<v|v>
# so the only thing ever computed with the state is the first column c_0 ... c_N, from ONE time
# evolution of |E> for N steps. No two states are ever added, so |E> can be an MPS as long as
# the time evolution is cheap (the backend is brickwall_tev through floquet_step).
#
# Notation: k = order of the geometric sum, ncv = number of Krylov vectors
# (Chapter 7 calls them m and k). The j-th Krylov vector has (j-1)k+1 coefficients.
###


#####################################
# The first column of the overlap matrix: one time evolution
#####################################

# <a|b> for the state types the time evolution can return (DenseSwitch turns an MPS into an ITensor)

inner_state(a::MPS, b::MPS)=inner(a, b)
inner_state(a::ITensor, b::ITensor)=(dag(a)*b)[]
inner_state(a::AbstractVector, b::AbstractVector)=dot(a, b)

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


#############
# Restoring the eigenvectors (as coefficients in the power basis) from the eigenvectors W of H
###########

restore_eigvecs(V::AbstractMatrix, W::AbstractMatrix)=V*W


#############
# The eigenvalues of U, from <psi|U|psi>/<psi|psi> with U|j> = |j+1>, i.e <i|U|j> = c_{j+1-i}
###########

function eigenvalues(Y::AbstractMatrix, columneinz)
    n=size(Y, 1)
    S=overlap(columneinz[1:n])
    S1=overlap_shift(columneinz, n)
    return [dot(y, S1*y)/real(dot(y, S*y)) for y in eachcol(Y)]
end


#############
# Convergence: ||U|psi> - omega|psi>|| for normalised |psi>
#
# Since U is unitary, ||U psi - omega psi||^2 = <psi|psi> (1 - |omega|^2) for omega = <psi|U|psi>/<psi|psi>,
# so the residual follows from the eigenvalue alone (and is resolved down to ~ sqrt(machine precision))
###########

# abs: in floating point a Rayleigh quotient can come out with |omega| > 1 (badly conditioned vectors),
# which is just as far from an eigenpair
convergence(omega)=sqrt.(abs.(1 .- abs2.(omega)))


#############
# Conditioning: kappa = ||y||^2 / <psi|psi>, coefficient norm against state norm
#
# The power basis is overcomplete (N > d) and its overlap matrix singular, so a vector can have
# coefficients much larger than its norm. Every quadratic form above is then computed with a
# relative error ~ kappa * (error of the c_j), which can fake |omega| ~ 1. An eigenpair is only
# trusted if kappa * overlap_error < tol as well.
###########

function conditioning(Y::AbstractMatrix, columneinz)
    n=size(Y, 1)
    S=overlap(columneinz[1:n])
    return [norm(y)^2/real(dot(y, S*y)) for y in eachcol(Y)]
end


#############
# POLFED in the power basis
###########

# U: FloquetOperator (or a matrix), psi0: |E> (MPS, dense ITensor or Vector)
# k: order of the geometric sum, ncv: number of Krylov vectors (k*ncv has to exceed the Hilbert space dimension)
# phi: target, tol: residual, overlap_error: accuracy of the c_j (machine precision for an exact time
# evolution, ~ sqrt of the discarded weight for a truncated one)
# returns omega, residuals and kappa (converged only, closest to the target first), the coefficients of
# the eigenvectors in the power basis, the overlaps c_j, and the number of Floquet steps (= ncv*k + 1)

#
# c: overlaps from an earlier run (length >= ncv*k+2) to skip the time evolution; one evolution of N
# steps serves every (k, ncv) with ncv*k+1 <= N

function polfed_mps(U, psi0; k::Int, ncv::Int, phi::Real=0.0, tol::Real=1e-5, overlap_error::Real=1e-14,
                    c::Union{Nothing,AbstractVector}=nothing)

    N=ncv*k+1
    if c===nothing
        c=overlaps(U, psi0, N)
    else
        length(c) >= N+1 || error("polfed_mps: need at least ncv*k+2 = $(N+1) overlaps, got $(length(c))")
        c=c[1:N+1]
    end

    H, V = arnoldi_power(c, k, ncv, phi)
    ritz, W = eigen(H)
    Y=restore_eigvecs(V, W)

    omega=eigenvalues(Y, c)
    res=convergence(omega)
    kap=conditioning(Y, c)

    conv=findall((res .< tol) .& (kap .* overlap_error .< tol))
    conv=conv[sortperm(abs.(omega[conv] .- exp(im*phi)))]

    return (omega=omega[conv], residuals=res[conv], kappa=kap[conv], coeffs=Y[:, conv], c=c, nsteps=N)
end


#############
# For checks: the eigenvector as a state, sum_j y_j U^j|E> (dense states only, one more time evolution)
###########

function power_state(U, psi0, y::AbstractVector)
    psi=psi0
    acc=y[1]*psi0
    for j in 2:length(y)
        psi=floquet_step(U, psi)
        acc=acc + y[j]*psi
    end
    return acc/norm(acc)
end
