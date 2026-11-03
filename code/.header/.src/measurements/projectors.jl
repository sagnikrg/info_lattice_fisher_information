# Z-basis projectors are diagonal: zero the amplitudes whose bit on site ind is the other value.
# reshape shares memory with Psi, so this is in place (no permutedims / matmul copies).
function zero_site!(Psi::Vector{ComplexF64}, ind::Int, s::Int)   # s=1 keeps ↑, s=2 keeps ↓

 L=Int(log2(length(Psi)))
 v=reshape(Psi, (2^(L-ind), 2, 2^(ind-1)))   # sites ind+1..L | site ind | sites 1..ind-1 (kron order: site 1 most significant)
 v[:, 3-s, :] .= 0

return Psi

end


# zero weight on this outcome: return the zero vector instead of 0/0 = NaN
function normalize_or_zero!(Psi::Vector{ComplexF64})

nrm=norm(Psi)
nrm>1e-12 ? rmul!(Psi, 1/nrm) : fill!(Psi, 0)

return Psi

end


function project_site!(Psi::Vector{ComplexF64}, ind::Int, s::Int)
 zero_site!(Psi, ind, s)
 return normalize_or_zero!(Psi)
end




function projector_down(Psi::Vector{ComplexF64}, ind::Int)
 return project_site!(copy(Psi), ind, 2)
end




function projector_up(Psi::Vector{ComplexF64}, ind::Int)
 return project_site!(copy(Psi), ind, 1)
end


function projectatt(Psi::Vector{ComplexF64}, indlist::Array{Int})

Psi=copy(Psi)
for n in indlist
    zero_site!(Psi, n, 2)
end

return normalize_or_zero!(Psi)

end

;
