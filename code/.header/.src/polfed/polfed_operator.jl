######################################################
# The Floquet unitary as an operator
#
# POLFED only ever needs U|psi>, i.e one Floquet step of the time evolution.
# The step is done by brickwall_tev (itensor/itev.jl), so the state decides the algorithm:
#   dense ITensor        : exact
#   MPS                  : with trunc = cutoff::Real, TruncationBudget or DenseSwitch
# For tests a plain matrix U acting on a Vector works as well.
######################################################


struct FloquetOperator{T}
    gatesEven::Vector{ITensor}
    gatesOdd::Vector{ITensor}
    trunc::T                    # nothing (exact), cutoff::Real, TruncationBudget or DenseSwitch
end

FloquetOperator(gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor})=FloquetOperator(gatesEven, gatesOdd, nothing)


# one Floquet step U|psi>

floquet_step(U::FloquetOperator{Nothing}, psi)=brickwall_tev(psi, U.gatesEven, U.gatesOdd)
floquet_step(U::FloquetOperator, psi)=brickwall_tev(psi, U.gatesEven, U.gatesOdd, U.trunc)

# dense matrix (small systems, tests)
floquet_step(U::AbstractMatrix, psi::AbstractVector)=U*psi


######################################################
# The geometric sum filter
#
#   g_k(U)|psi> = sum_{m=0}^{k} e^{-i m phi} U^m |psi>
#
# k Floquet steps, summed up with the phases. The sum of states only stays cheap for dense
# states (Vector or dense ITensor); for an MPS the bond dimensions add up, which is why
# polfed_mps.jl never sums states.
######################################################

# brickwall_tev renormalises after every step, so the state is evolved normalised and the
# norm is put back at the end (keeps g_k linear for any input, e.g the vectors ARPACK sends)

function geometric_sum(U, psi, k::Int, phi::Real)
    nrm=norm(psi)
    cur=psi/nrm
    acc=copy(cur)
    for m in 1:k
        cur=floquet_step(U, cur)
        acc=acc + exp(-im*m*phi)*cur
    end
    return nrm*acc
end

# the same on the unit circle, g_k(e^{i phi_n}) (Luitz Eq. (8)), for plots and checks
geometric_sum(z::Number, k::Int, phi::Real)=sum(exp(-im*m*phi)*z^m for m in 0:k)


######################################################
# The full Floquet matrix, column by column from the basis states (exact, dense evolution)
# for full diagonalisation as the reference (small L), in the index order of sites
######################################################

function floquet_matrix(U::FloquetOperator, sites::Vector{<:Index})
    L=length(sites)
    d=2^L
    M=zeros(ComplexF64, d, d)
    for i in 1:d
        e=zeros(ComplexF64, d)
        e[i]=1.0
        psi=ITensor(reshape(e, fill(2, L)...), sites...)
        M[:, i]=vec(Array(brickwall_tev(psi, U.gatesEven, U.gatesOdd), sites...))
    end
    return M
end
