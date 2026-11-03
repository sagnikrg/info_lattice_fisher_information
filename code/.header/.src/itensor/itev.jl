

# gates typed so these do not overwrite brickwall_tev(Psi, brick, sites, dummysites) in circuit/brickwall_tev.jl

function brickwall_tev(psi, gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor})

  psi = apply(gatesEven, psi)
  psi = apply(gatesOdd, psi)
  normalize!(psi)           # Normalize the state after applying gates

  return psi
end

# dense ITensor: no SVDs to truncate, cutoff is ignored (apply on ITensor does not accept it)
function brickwall_tev(psi, gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor}, cutoff::Real)

  brickwall_tev(psi, gatesEven, gatesOdd)
end

# MPS: truncated evolution
function brickwall_tev(psi::MPS, gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor}, cutoff::Real)

  psi = apply(gatesEven, psi; cutoff)
  psi = apply(gatesOdd, psi; cutoff)
  normalize!(psi)           # Normalize the state after applying gates

  return psi
end


#-------------------------------------
# dynamic truncation 1: error budget
#-------------------------------------
#
# Instead of a fixed cutoff, the whole run gets a budget delta of discarded weight
# (sum of the discarded squared Schmidt values over all SVDs, ~ the infidelity of the final state).
# Before every step the cutoff is set to (budget left)/(SVDs left), and after the step the weight
# actually discarded is read off from the norm loss of the MPS (apply does not renormalise,
# and every SVD is done at the orthogonality centre, so 1 - <psi|psi> is exactly what was cut).
# Steps that need to discard less than allowed (small chi, early times) save budget for later.
#
#   budget=TruncationBudget(delta, tnum, L)
#   for t in 1:tnum; psi=brickwall_tev(psi, gatesEven, gatesOdd, budget); end

mutable struct TruncationBudget
    delta::Float64          # total discarded weight allowed over the run
    nsvd::Int               # SVDs left in the run
    spent::Float64          # discarded weight so far
    cutoffs::Vector{Float64}    # cutoff used in every step (for diagnostics)
end

TruncationBudget(delta::Real, tnum::Int, L::Int)=TruncationBudget(delta, tnum*(L-1), 0.0, Float64[])

# dense ITensor: exact, nothing to budget
function brickwall_tev(psi, gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor}, budget::TruncationBudget)

  brickwall_tev(psi, gatesEven, gatesOdd)
end

function brickwall_tev(psi::MPS, gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor}, budget::TruncationBudget)

  nstep=length(gatesEven)+length(gatesOdd)
  cutoff=max(budget.delta-budget.spent, 0.0)/max(budget.nsvd, nstep)
  push!(budget.cutoffs, cutoff)

  psi = apply(gatesEven, psi; cutoff)
  psi = apply(gatesOdd, psi; cutoff)

  budget.spent+=max(1-real(inner(psi, psi)), 0.0)    # weight discarded in this step
  budget.nsvd-=nstep

  normalize!(psi)           # Normalize the state after applying gates

  return psi
end


#-------------------------------------
# dynamic truncation 2: switch to dense
#-------------------------------------
#
# Evolve as an MPS (with cutoff) while it is cheap, and contract it into one dense L legged
# ITensor once the bond dimension reaches chimax, from where on the evolution is exact.
# The return type changes from MPS to ITensor at the switch, multiple dispatch takes care of the rest.
#
#   switch=DenseSwitch(cutoff, chimax)
#   for t in 1:tnum; psi=brickwall_tev(psi, gatesEven, gatesOdd, switch); end

struct DenseSwitch
    cutoff::Float64         # cutoff while still an MPS
    chimax::Int             # switch to dense once maxlinkdim(psi) >= chimax
end

# already dense: exact
function brickwall_tev(psi, gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor}, switch::DenseSwitch)

  brickwall_tev(psi, gatesEven, gatesOdd)
end

function brickwall_tev(psi::MPS, gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor}, switch::DenseSwitch)

  psi = brickwall_tev(psi, gatesEven, gatesOdd, switch.cutoff)

  maxlinkdim(psi) >= switch.chimax ? prod(psi) : psi
end
