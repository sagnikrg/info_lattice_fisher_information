
include("idtcckt.jl")


function Brickwall_tev(psi,gatesEven,gatesOdd; method=:ITensor)

  psi = apply(gatesEven, psi)
  psi = apply(gatesOdd, psi)
  normalize!(psi)           # Normalize the state after applying gates

  return psi
end

function Brickwall_tev(psi,gatesEven,gatesOdd, cutoff; method=:ITensor)

  psi = apply(gatesEven, psi; cutoff)
  psi = apply(gatesOdd, psi; cutoff)
  normalize!(psi)           # Normalize the state after applying gates

  return psi
end