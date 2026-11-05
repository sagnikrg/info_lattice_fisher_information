### Sector resolved eigendecomposition of a unitary
###
### Every function here returns nested arrays, one entry per sector:
###
###     E, V = eig_ztot(U)
###     E[k]        eigenvalues of sector k
###     V[k]        eigenvectors of sector k as columns, in the full 2^L basis
###
### The order of the sectors is that of symmetry_sectors(L; ...), which gives the labels.
### vectors=false skips the eigenvectors, V is then nothing.

using LinearAlgebra


# Schur instead of eigen: U is normal, so T is diagonal and the vectors stay orthonormal inside degeneracies

function eig_block(B, vectors)

    if vectors
        F=schur(B)
        return F.values, F.vectors
    end

    return eigvals(B), nothing
end


"""
    eig_sectors(U, Qs; vectors=true, tol=1e-8)

Blocks `Q'*U*Q` for every isometry in `Qs`, and the eig of each block.
Throws an error if `U` leaks out of a sector, i.e. if a block is not unitary.
"""
function eig_sectors(U, Qs; vectors=true, tol=1e-8)

    E=Vector{Vector{ComplexF64}}(undef, length(Qs))
    V=vectors ? Vector{Matrix{ComplexF64}}(undef, length(Qs)) : nothing

    for (k,Q) in enumerate(Qs)

        B=Matrix{ComplexF64}(Q'*U*Q)

        leak=norm(B'*B-I)/sqrt(size(B,1))
        leak<tol || error("U is not block diagonal in sector $k: |B'B-1|/sqrt(d) = $leak")

        values, vecs = eig_block(B, vectors)

        E[k]=values
        vectors && (V[k]=Q*vecs)
    end

    return E, V
end


"""
    eig_symmetry(U; ztot, sz, spin_parity, chain_parity, highest_weight=false, vectors=true, tol=1e-8)

`eig_sectors` for any combination of the symmetries in `symmetry_sectors`.
"""
function eig_symmetry(U; vectors=true, tol=1e-8, symmetries...)

    L=Int(log2(size(U,1)))
    Qs, _ = symmetry_sectors(L; symmetries...)

    return eig_sectors(U, Qs; vectors=vectors, tol=tol)
end


# Single symmetry wrappers

eig_ztot(U; kwargs...)         = eig_symmetry(U; ztot=true, kwargs...)
eig_sz(U; kwargs...)           = eig_symmetry(U; sz=true, kwargs...)
eig_spin_parity(U; kwargs...)  = eig_symmetry(U; spin_parity=true, kwargs...)
eig_chain_parity(U; kwargs...) = eig_symmetry(U; chain_parity=true, kwargs...)
eig_abs_ztot(U; kwargs...)     = eig_symmetry(U; abs_ztot=true, kwargs...)


"""
    eig_perfect_kick(U; tol=1e-8)

For a perfect kick, U maps m to -m, so that levels of |m|>0 come in exact pairs (E,-E) and are two copies of one spectrum.
Returns one nested list of phases per |m|, each spectrum once: eigenvalues of U^2 inside m for m>0 and of U inside m=0.
For odd L there is no m=0. Only for statistics of the levels, the eigenvalues of U^2 are 2*phi of the pairs.
"""
function eig_perfect_kick(U; tol=1e-8)

    L=Int(log2(size(U,1)))

    Qs, labels = symmetry_sectors(L; ztot=true)
    keep=[l.m>0 for l in labels]
    E, _ = eig_sectors(U^2, Qs[keep]; vectors=false, tol=tol)

    if iseven(L)
        Q0=Qs[findfirst(l->l.m==0, labels)]
        E0, _ = eig_sectors(U, [Q0]; vectors=false, tol=tol)
        append!(E, E0)
    end

    return E
end
