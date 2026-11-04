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

# the building blocks are in functions/power_basis.jl (overlaps, inner_prod, apply_polynomial, arnoldi_power)
# and functions/power_eigenpairs.jl (restore_eigvecs, eigenvalues, convergence, conditioning)


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
# The pi pairs of the time crystal: eigenvalues around phi and around phi + pi from the same overlaps
#
# The target only enters the Krylov part (the phases of g_k), so one time evolution serves any target.
# In the time crystal the eigenstates come in pairs with quasienergies epsilon and epsilon + pi
# (+ an exponentially small splitting), so the second target gives the partners of the first. Not
# possible in the normal POLFED without a second run, since there the Krylov vectors are states.
# returns both runs (near: around phi, far: around phi + pi), and for every eigenvalue omega of near
# its partner (the eigenvalue of far closest to -omega) and the pi pair splitting |arg(-partner/omega)|
###########

function polfed_mps_pairs(U, psi0; k::Int, ncv::Int, phi::Real=0.0, c::Union{Nothing,AbstractVector}=nothing, kwargs...)

    near=polfed_mps(U, psi0; k=k, ncv=ncv, phi=phi, c=c, kwargs...)
    far=polfed_mps(U, psi0; k=k, ncv=ncv, phi=phi+pi, c=near.c, kwargs...)

    partner=ComplexF64[]; splitting=Float64[]
    for w in near.omega
        isempty(far.omega) && break
        p=far.omega[argmin(abs.(far.omega .+ w))]
        push!(partner, p); push!(splitting, abs(angle(-p/w)))
    end

    return (near=near, far=far, partner=partner, splitting=splitting)
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
