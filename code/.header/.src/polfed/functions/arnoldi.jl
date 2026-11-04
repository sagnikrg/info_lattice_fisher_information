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
