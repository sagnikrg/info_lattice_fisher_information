######################################################
# Normal POLFED: Arnoldi on g_k(U) with explicit Krylov vectors
#
# The Krylov vectors are dense states (Vector or dense L legged ITensor), every application of
# g_k(U) is k Floquet steps (geometric_sum, polfed_operator.jl).
#
#   method = :arpack  : implicitly restarted Arnoldi (ARPACK, as in Luitz 2021), g_k(U) as a matrix free operator
#   method = :arnoldi : plain Arnoldi with ncv vectors, no restarts (simple, needs a larger ncv)
#
# The eigenvalues of U are read off from the eigenvectors, omega_n = <n|U|n>, with the residual
# r_n = ||U|n> - omega_n|n>|| as the convergence check (Luitz Eq. (10)).
######################################################


######################################################
# State <-> Vector (ARPACK works on plain vectors)
######################################################

# the index order is fixed by the template (the time evolution may permute the indices of a dense ITensor)

to_vector(psi::AbstractVector, template::AbstractVector=psi)=psi
to_vector(psi::ITensor, template::ITensor=psi)=vec(Array(psi, inds(template)...))

from_vector(v::AbstractVector, template::AbstractVector)=v
from_vector(v::AbstractVector, template::ITensor)=ITensor(reshape(v, dim.(inds(template))...), inds(template)...)


######################################################
# Arnoldi on a linear map f with full re-orthogonalisation (classical Gram-Schmidt, twice)
# works for any state with dot, norm, + and scalar * (Vector, dense ITensor)
######################################################

function arnoldi(f::Function, v0, ncv::Int)

    V=Any[v0/norm(v0)]
    H=zeros(ComplexF64, ncv+1, ncv)

    for j in 1:ncv
        w=f(V[j])
        for pass in 1:2, i in 1:j                  # second pass restores the orthogonality lost to round off
            h=dot(V[i], w)
            H[i, j]+=h
            w=w - h*V[i]
        end
        H[j+1, j]=norm(w)
        if abs(H[j+1, j]) < 1e-14                  # invariant subspace found
            return H[1:j, 1:j], V[1:j]
        end
        j<ncv && push!(V, w/H[j+1, j])
    end

    return H[1:ncv, 1:ncv], V
end


######################################################
# Eigenpairs of U from approximate eigenvectors
######################################################

# omega = <x|U|x> and r = ||U|x> - omega|x>|| for normalised x, sorted by distance to the target

function unitary_eigenpairs(U, X, phi::Real; tol::Real=1e-8)

    X=[x/norm(x) for x in X]
    UX=[floquet_step(U, x) for x in X]
    omega=[dot(x, ux) for (x, ux) in zip(X, UX)]
    res=[norm(ux - w*x) for (x, ux, w) in zip(X, UX, omega)]

    conv=findall(res .< tol)
    conv=conv[sortperm(abs.(omega[conv] .- exp(im*phi)))]

    return (omega=omega[conv], vectors=X[conv], residuals=res[conv])
end


######################################################
# g_k(U) as a matrix free operator for ARPACK
######################################################

mutable struct GeometricSumMap{TU,TS} <: AbstractMatrix{ComplexF64}
    U::TU
    template::TS                # a state, to rebuild states from vectors (dense ITensor or Vector)
    k::Int
    phi::Float64
    nmatvec::Int                # applications of g_k(U) so far
end

Base.size(G::GeometricSumMap)=(d=length(to_vector(G.template)); (d, d))
Base.getindex(G::GeometricSumMap, i::Int, j::Int)=error("GeometricSumMap is matrix free")

function LinearAlgebra.mul!(y::AbstractVector, G::GeometricSumMap, x::AbstractVector)
    G.nmatvec+=1
    y .= to_vector(geometric_sum(G.U, from_vector(Vector{ComplexF64}(x), G.template), G.k, G.phi), G.template)
    return y
end


######################################################
# POLFED
######################################################

# U    : FloquetOperator (or a matrix), psi0: initial dense state (Vector or dense ITensor)
# k    : order of the geometric sum, ncv: Krylov space dimension, nev: requested eigenpairs (arpack)
# phi  : target phase, tol: residual for an eigenpair to count as converged
# returns omega, vectors, residuals (converged only, closest to the target first),
# and the number of applications of g_k(U) and of Floquet steps

function polfed(U, psi0; k::Int, ncv::Int, nev::Int=div(ncv, 2), phi::Real=0.0, tol::Real=1e-8,
                method::Symbol=:arpack, arpack_tol::Real=1e-12, maxiter::Int=300)

    if method==:arnoldi

        nmatvec=Ref(0)
        f(v)=(nmatvec[]+=1; geometric_sum(U, v, k, phi))
        H, V = arnoldi(f, psi0, ncv)
        ritz, W = eigen(H)
        X=[sum(W[i, n]*V[i] for i in eachindex(V)) for n in axes(W, 2)]
        nmv=nmatvec[]

    elseif method==:arpack

        G=GeometricSumMap(U, psi0, k, Float64(phi), 0)
        ritz, Xm = Arpack.eigs(G; nev=nev, ncv=ncv, which=:LM, tol=arpack_tol, maxiter=maxiter, v0=Vector{ComplexF64}(to_vector(psi0)))
        X=[from_vector(Xm[:, n], psi0) for n in axes(Xm, 2)]
        nmv=G.nmatvec

    else
        error("polfed: method should be :arpack or :arnoldi")
    end

    pairs=unitary_eigenpairs(U, X, phi; tol=tol)

    return (omega=pairs.omega, vectors=pairs.vectors, residuals=pairs.residuals,
            nmatvec=nmv, nsteps=nmv*k + length(X))
end
