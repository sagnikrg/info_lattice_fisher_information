######################################################
# Header file dedicated to polynomial filtering
#
# POLFED: Polynomial Filtering Exact Diagonalisation (Luitz, SciPost Phys. 11, 021 (2021))
# Eigenpairs of the Floquet unitary U closest to a target e^{i phi} on the unit circle,
# from the Arnoldi algorithm on the geometric sum g_k(U) = sum_{m=0}^{k} e^{-i m phi} U^m.
#
#   polfed.jl     : normal POLFED (this file), Arnoldi on g_k(U) with explicit Krylov vectors
#   polfed_mps.jl : POLFED in the power basis U^j|E>, only the overlaps <E|U^j|E> are needed,
#                   i.e one time evolution, so |E> can be an MPS (Chapter 7)
#
#   functions/operator.jl         : FloquetOperator (itev backend), floquet_step, geometric_sum, floquet_matrix
#   functions/states.jl           : to_vector, from_vector, inner_state
#   functions/arnoldi.jl          : arnoldi (linear map), unitary_eigenpairs, GeometricSumMap (for ARPACK)
#   functions/overlap.jl          : overlap, overlap_shift (Toeplitz overlap matrices of the power basis)
#   functions/power_basis.jl      : overlaps, inner_prod, apply_polynomial, arnoldi_power
#   functions/power_eigenpairs.jl : restore_eigvecs, eigenvalues, convergence, conditioning
#   legacy/                       : the first matrix version (not included here)
#
# Both use brickwall_tev as the time evolution backend through FloquetOperator: normal POLFED needs dense
# Krylov vectors (Vector or dense ITensor), the power basis takes any state (dense, MPS, truncated, switched).
# needs itensor/ITensor.jl (brickwall_tev, TruncationBudget, DenseSwitch) included before
######################################################

using LinearAlgebra
using Arpack

include("functions/operator.jl")
include("functions/states.jl")
include("functions/arnoldi.jl")
include("functions/overlap.jl")
include("functions/power_basis.jl")
include("functions/power_eigenpairs.jl")


######################################################
# Normal POLFED: Arnoldi on g_k(U) with explicit Krylov vectors
#
# The Krylov vectors are dense states (Vector or dense L legged ITensor), every application of g_k(U)
# is k Floquet steps (geometric_sum, functions/operator.jl).
#
#   method = :arpack  : implicitly restarted Arnoldi (ARPACK, as in Luitz 2021), g_k(U) as a matrix free operator
#   method = :arnoldi : plain Arnoldi with ncv vectors, no restarts (simple, needs a larger ncv)
#
# The eigenvalues of U are read off from the eigenvectors, omega_n = <n|U|n>, with the residual
# r_n = ||U|n> - omega_n|n>|| as the convergence check (Luitz Eq. (10)).
######################################################


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


include("polfed_mps.jl")
