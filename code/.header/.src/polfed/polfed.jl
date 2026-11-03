######################################################
# Header file dedicated to polynomial filtering
#
# POLFED: Polynomial Filtering Exact Diagonalisation (Luitz, SciPost Phys. 11, 021 (2021))
# Eigenpairs of the Floquet unitary U closest to a target e^{i phi} on the unit circle,
# from the Arnoldi algorithm on the geometric sum g_k(U) = sum_{m=0}^{k} e^{-i m phi} U^m.
#
#   polfed_operator.jl : the Floquet unitary as an operator (FloquetOperator, itev backend) and g_k(U)|psi>
#   polfed_dense.jl    : normal POLFED, Arnoldi on g_k(U) with explicit Krylov vectors (dense states)
#   polfed_mps.jl      : POLFED in the power basis U^j|E>, only the overlaps <E|U^j|E> are needed,
#                        i.e one time evolution, so |E> can be an MPS (Chapter 7)
#
# needs itensor/ITensor.jl (brickwall_tev, TruncationBudget, DenseSwitch) included before
######################################################

using LinearAlgebra
using Arpack

include("functions/overlap.jl")
include("functions/mat_mult.jl")


######################################################
# Base code with matrix multiplication
######################################################


# Function to perform Arnoldi iteration
function arnoldi(A::AbstractMatrix, v::AbstractVector, k::Int)
    n = length(v)
    H = zeros(ComplexF64, k, k)  # Hessenberg matrix
    V = zeros(ComplexF64, n, k)  # Orthonormal basis

    V[:, 1] = v / norm(v)  # Normalize the initial vector

    for j in 1:k                     # all k columns, so that H = V' A V is complete (exact for k = n)
        w = A * V[:, j]  # Apply matrix A to the j-th basis vector

        for i in 1:j
            H[i, j] = dot(V[:, i], w)  # Compute the j-th column of H
            w -= H[i, j] * V[:, i]  # Orthogonalize w against previous basis vectors
        end

        j == k && break              # H[k+1, k] would lie outside the k x k matrix
        H[j + 1, j] = norm(w)  # Compute the (j+1)-th entry of the j-th column of H

        if H[j + 1, j] != 0
            V[:, j + 1] = w / H[j + 1, j]  # Normalize the (j+1)-th basis vector
        end
    end

    return H, V
end

# Function to perform Arnoldi diagonalization
function arnoldi_diagonalize(A::AbstractMatrix, k::Int)
    n = size(A, 1)
    v0 = rand(ComplexF64, n)  # Random initial vector
    H, V = arnoldi(A, v0, k)  # Perform Arnoldi iteration

    # Diagonalize H
    evals, evecs = eigen(H[1:k, 1:k])

    # Approximate eigenvectors of A
    approx_eigenvectors = V[:, 1:k] * evecs

    return evals, approx_eigenvectors
end



# Rayleigh quotients <v|A|v>/<v|v> of the columns of V

function extract_eigenvalues(A, V)
    k = size(V, 2)
    eigenvalues = zeros(ComplexF64, k)
    for i in 1:k
        v = V[:, i]
        Av = A * v
        eigenvalue = dot(v, Av) / dot(v, v)  # Rayleigh quotient
        eigenvalues[i] = eigenvalue
    end
    return eigenvalues
end


# And finally to sort the best ones: the n largest in modulus (default: half of them)

function filter_eig(eigvalsfull, n::Int=div(length(eigvalsfull), 2))
    eigvalsfull=sort(eigvalsfull, by=abs, rev=true)
    eigvalsfull = eigvalsfull[1:n]
    return eigvalsfull
end


include("polfed_operator.jl")
include("polfed_dense.jl")
include("polfed_mps.jl")
