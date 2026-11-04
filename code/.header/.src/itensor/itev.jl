

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


#-------------------------------------
# adaptive time evolution: the cheapest representation at every step
#-------------------------------------
#
# The default for production runs. The state is evolved as an MPS (with cutoff) as long as an MPS
# step is cheaper than a dense one, and as a dense L legged ITensor (exact) from then on:
#   - the cost of a dense step at size L is measured once, on the first call (from prod(psi), while
#     the MPS is still cheap to contract), and kept in the struct, so one AdaptiveTEV serves all runs
#     of the same L in a job
#   - every MPS step is timed, and once nslow steps in a row were slower than a dense step the state
#     is contracted (the bond dimension only grows, so the switch is one way)
#   - cores: the run may use up to ncores cores (BLAS and Strided threads). They are only switched on
#     while the MPS bond dimension is >= chithreads, where the SVDs are large enough to gain (~2x on 4
#     cores), everything else runs on one thread (the dense step does not gain from threads in ITensor)
# Accuracy: the infidelity is bounded by the cutoff, t(L-1)*cutoff for the t MPS steps, and exact after.
# The switching time depends on the machine (it is timed), the result does not beyond that bound.
#
#   tev=AdaptiveTEV(; cutoff=1e-10, ncores=1)
#   for run in runs
#       reset!(tev); psi=random_bitstring(sites)                  # an MPS
#       for t in 1:tnum; psi=brickwall_tev(psi, gatesEven, gatesOdd, tev); end
#   end
#   tev.history, tev.chi : representation (:mps, :dense) and bond dimension after every step

using LinearAlgebra: BLAS
using ITensors.NDTensors: Strided

mutable struct AdaptiveTEV
    cutoff::Float64             # cutoff while an MPS
    ncores::Int                 # cores the run may use
    chithreads::Int             # use ncores threads from this bond dimension on
    nslow::Int                  # switch after this many MPS steps in a row slower than a dense step
    tdense::Dict{Int,Float64}   # L => wall time of one dense step (single thread)
    slow::Int                   # MPS steps in a row slower than dense so far
    history::Vector{Symbol}     # :mps or :dense, for every step
    chi::Vector{Int}            # bond dimension after every step (0 when dense)
end

AdaptiveTEV(; cutoff::Real=1e-10, ncores::Int=1, chithreads::Int=500, nslow::Int=2)=
    AdaptiveTEV(cutoff, ncores, chithreads, nslow, Dict{Int,Float64}(), 0, Symbol[], Int[])

# a new run with the same calibration
function reset!(tev::AdaptiveTEV)
    tev.slow=0
    empty!(tev.history)
    empty!(tev.chi)
    return tev
end

function set_cores!(n::Int)
    BLAS.set_num_threads(n)
    Strided.set_num_threads(min(n, Threads.nthreads()))
end

# time of one dense step at the size of psi, single thread, fastest of two (the first may compile)
function calibrate!(tev::AdaptiveTEV, psi::MPS, gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor})
    set_cores!(1)
    phi=prod(psi)
    t=minimum(@elapsed(brickwall_tev(phi, gatesEven, gatesOdd)) for _ in 1:2)
    tev.tdense[length(psi)]=t
    return t
end

# already dense: exact, single thread
function brickwall_tev(psi, gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor}, tev::AdaptiveTEV)

  set_cores!(1)
  psi = brickwall_tev(psi, gatesEven, gatesOdd)
  push!(tev.history, :dense); push!(tev.chi, 0)

  return psi
end

function brickwall_tev(psi::MPS, gatesEven::Vector{ITensor}, gatesOdd::Vector{ITensor}, tev::AdaptiveTEV)

  L=length(psi)
  haskey(tev.tdense, L) || calibrate!(tev, psi, gatesEven, gatesOdd)

  set_cores!(maxlinkdim(psi) >= tev.chithreads ? tev.ncores : 1)
  t = @elapsed psi = brickwall_tev(psi, gatesEven, gatesOdd, tev.cutoff)
  push!(tev.history, :mps); push!(tev.chi, maxlinkdim(psi))

  tev.slow = t > tev.tdense[L] ? tev.slow + 1 : 0
  if tev.slow >= tev.nslow
    set_cores!(1)
    return prod(psi)
  end

  return psi
end
