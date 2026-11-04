#######################################
# Eigenpairs of U from the Ritz vectors in the power basis: eigenvectors, eigenvalues, residuals, conditioning
# (used by polfed_mps.jl)
#######################################


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
