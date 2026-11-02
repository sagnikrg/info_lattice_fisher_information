

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
