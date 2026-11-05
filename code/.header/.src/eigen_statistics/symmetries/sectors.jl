### Bases of the symmetry sectors of a spin-1/2 chain of length L
###
### ztot         : total magnetisation m, eigenvalue of sum_j Z_j / 2
### sz           : total spin and magnetisation |s,m>, eigenvalues of S^2 and sum_j Z_j / 2
### spin_parity  : p_z = +-1, eigenvalue of prod_j Z_j
### chain_parity : p_c = +-1, eigenvalue of the reflection R, site j -> L+1-j
###
### All of them commute with each other, so they can be combined freely.
### A sector is stored as an isometry Q (2^L x d_sector), block = Q'*U*Q.

using LinearAlgebra
using SparseArrays


################################
# Operators
################################


# Magnetisation of the basis state with index i (1 based): up=|0>, down=|1>

magnetisation(i, L) = (L - 2*count_ones(i-1))/2

# Eigenvalue of prod_j Z_j on the basis state with index i

spin_parity_value(i) = iseven(count_ones(i-1)) ? 1 : -1

# Basis state with index i mirrored, site j -> L+1-j

reflect_index(i, L) = 1 + Int(bitreverse(UInt64(i-1)) >> (64-L))


# S^2 = 3L/4 - L(L-1)/4 + sum_{i<j} P_ij , with P_ij the swap of sites i and j

function total_spin_squared(L)

    dim=2^L
    rows=Int[]; cols=Int[]; vals=Float64[]
    diagonal=fill(3L/4 - L*(L-1)/4, dim)

    for k in 0:dim-1, i in 0:L-2, j in i+1:L-1
        if (k>>i)&1 == (k>>j)&1
            diagonal[k+1]+=1                     # aligned, swap does nothing
        else
            push!(rows, k+1)
            push!(cols, xor(k, (1<<i)|(1<<j))+1)
            push!(vals, 1.0)
        end
    end

    return sparse(rows, cols, vals, dim, dim) + spdiagm(0 => diagonal)
end



################################
# Sectors
################################


# Basis states grouped by the quantum numbers that are diagonal in the Z basis

function diagonal_sectors(L, resolve_m, resolve_pz, resolve_absm=false)

    groups=Dict{NamedTuple, Vector{Int}}()

    for i in 1:2^L
        label=NamedTuple()
        resolve_m  && (label=merge(label, (m=magnetisation(i,L),)))
        resolve_absm && (label=merge(label, (absm=abs(magnetisation(i,L)),)))
        resolve_pz && (label=merge(label, (pz=spin_parity_value(i),)))
        push!(get!(groups, label, Int[]), i)
    end

    labels=sort(collect(keys(groups)), by=l->(-get(l,:m,0.0), -get(l,:absm,0.0), -get(l,:pz,0)))

    return [groups[l] for l in labels], labels
end


# Symmetric and antisymmetric combinations of mirror images inside one group of basis states

function chain_parity_bases(indices, L)

    dim=2^L
    rows_p=Int[]; cols_p=Int[]; vals_p=Float64[]
    rows_m=Int[]; cols_m=Int[]; vals_m=Float64[]
    np=0; nm=0

    for i in indices
        r=reflect_index(i,L)
        if r==i
            np+=1
            append!(rows_p, i); append!(cols_p, np); append!(vals_p, 1.0)
        elseif r>i
            np+=1; nm+=1
            append!(rows_p, [i,r]); append!(cols_p, [np,np]); append!(vals_p, [1.0, 1.0]/sqrt(2))
            append!(rows_m, [i,r]); append!(cols_m, [nm,nm]); append!(vals_m, [1.0,-1.0]/sqrt(2))
        end
    end

    return sparse(rows_p, cols_p, vals_p, dim, np), sparse(rows_m, cols_m, vals_m, dim, nm)
end


# Splits the sector Q by total spin s, from the highest s to the lowest

function split_total_spin(Q, S2)

    F=eigen(Hermitian(Matrix(Q'*S2*Q)))
    s=(sqrt.(max.(1 .+ 4*F.values, 0)) .- 1)/2

    maximum(abs.(2*s .- round.(2*s)))<1e-8 || error("S^2 eigenvalues are not of the form s(s+1)")

    # grouped as the integer 2s, a float s can come out as -0.0 and unique would count it apart from 0.0
    twice_s=round.(Int, 2*s)
    spins=sort(unique(twice_s), rev=true)

    return [sparse(Q*F.vectors[:, twice_s .== t]) for t in spins], spins/2
end



"""
    symmetry_sectors(L; ztot, abs_ztot, sz, spin_parity, chain_parity, highest_weight=false)

Isometries and labels of the sectors of the chosen symmetries, returns `(Qs, labels)`.
`sz=true` resolves |s,m> and already contains `ztot` and `spin_parity`.
`abs_ztot=true` resolves |m|, for U that maps m to -m (a perfect kick).
`highest_weight=true` keeps only s=m, one copy of every SU(2) multiplet.
"""
function symmetry_sectors(L::Int; ztot=false, abs_ztot=false, sz=false, spin_parity=false, chain_parity=false, highest_weight=false)

    highest_weight && !sz && error("highest_weight needs sz=true")

    # spin_parity is a function of m, so adding it to ztot or sz only adds a label
    groups, labels = diagonal_sectors(L, ztot||sz, spin_parity, abs_ztot)

    if highest_weight
        keep=[l.m>=0 for l in labels]
        groups, labels = groups[keep], labels[keep]
    end

    Qs=[sparse(g, 1:length(g), ones(length(g)), 2^L, length(g)) for g in groups]


    if chain_parity
        Qs_new=SparseMatrixCSC{Float64,Int}[]
        labels_new=NamedTuple[]
        for (g,l) in zip(groups, labels)
            plus, minus = chain_parity_bases(g, L)
            size(plus,2)>0  && (push!(Qs_new, plus);  push!(labels_new, merge(l, (pc=1,))))
            size(minus,2)>0 && (push!(Qs_new, minus); push!(labels_new, merge(l, (pc=-1,))))
        end
        Qs, labels = Qs_new, labels_new
    end


    if sz
        S2=total_spin_squared(L)
        Qs_new=SparseMatrixCSC{Float64,Int}[]
        labels_new=NamedTuple[]
        for (Q,l) in zip(Qs, labels)
            Q_spin, spins = split_total_spin(Q, S2)
            for (Qs_s, s) in zip(Q_spin, spins)
                (highest_weight && s != l.m) && continue
                push!(Qs_new, Qs_s)
                push!(labels_new, merge((s=s,), l))
            end
        end
        Qs, labels = Qs_new, labels_new
    end

    return Qs, labels
end
