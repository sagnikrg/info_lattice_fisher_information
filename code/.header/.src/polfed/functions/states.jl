######################################################
# State <-> Vector (ARPACK works on plain vectors)
######################################################

# the index order is fixed by the template (the time evolution may permute the indices of a dense ITensor)

to_vector(psi::AbstractVector, template::AbstractVector=psi)=psi
to_vector(psi::ITensor, template::ITensor=psi)=vec(Array(psi, inds(template)...))

from_vector(v::AbstractVector, template::AbstractVector)=v
from_vector(v::AbstractVector, template::ITensor)=ITensor(reshape(v, dim.(inds(template))...), inds(template)...)


######################################################
# Inner product <a|b> for the state types the time evolution can return
######################################################

# (DenseSwitch turns an MPS into an ITensor)

inner_state(a::MPS, b::MPS)=inner(a, b)
inner_state(a::ITensor, b::ITensor)=(dag(a)*b)[]
inner_state(a::AbstractVector, b::AbstractVector)=dot(a, b)
