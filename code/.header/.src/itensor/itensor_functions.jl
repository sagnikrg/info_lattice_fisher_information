
#-------------------------------------
# half chain EntanglementEntropy for itensor
#-------------------------------------

function halfchainee(Psi, sites::Vector{<:Index})

    L=length(sites)
    Lhalf=Int(L/2)
    N=2^Lhalf;


    lft=combiner(sites[1:Lhalf], tags="left half")
    rgt=combiner(sites[Lhalf+1:L], tags="right half")
    Psi_svd=Psi*lft
    Psi_svd=Psi_svd*rgt


    U,V,D=svd(Psi_svd, inds(Psi_svd)[1]);
    ee=0.0

for i in 1:N
    if V[i,i]>1e-10
    ee=ee+-V[i,i]^2*log(V[i,i]^2)
    end
end

ee
end

# MPS: Schmidt values from the bond between sites Lhalf and Lhalf+1

function halfchainee(psi::MPS, sites::Vector{<:Index})

    L=length(sites)
    Lhalf=Int(L/2)

    psi=orthogonalize(psi, Lhalf)
    U,S,V=svd(psi[Lhalf], uniqueinds(psi[Lhalf], psi[Lhalf+1]));
    ee=0.0

for i in 1:dim(S, 1)
    if S[i,i]>1e-10
    ee=ee+-S[i,i]^2*log(S[i,i]^2)
    end
end

ee
end

#-------------------------------------
# <Z_i(0) Z_j(t)> for itensor
#-------------------------------------

# <Z_i(0) Z_j(t)> = <psi_0| Z_i U^{†t} Z_j U^t |psi_0> = <phi_t| Z_j |psi_t>   (Pauli Z)
#
#   psi_t = U^t psi_0
#   phi_t = U^t Z_i psi_0      (prepare Z_i psi_0 with apply_Z, then evolve it like psi_t)
#
# complex in general

#----- helper: Z_i|psi>

# dense ITensor
function apply_Z(psi, i::Int, sites::Vector{<:Index})

    noprime(op("Z", sites[i])*psi)
end

# MPS
function apply_Z(psi::MPS, i::Int, sites::Vector{<:Index})

    apply(op("Z", sites[i]), psi)
end


#----- general case

# dense ITensor
function ZZ(phi_t, psi_t, j::Int, sites::Vector{<:Index})

    (dag(phi_t)*apply_Z(psi_t, j, sites))[]
end

# MPS
function ZZ(phi_t::MPS, psi_t::MPS, j::Int, sites::Vector{<:Index})

    inner(phi_t, apply_Z(psi_t, j, sites))
end

# same site: <Z_i(0) Z_i(t)>, phi_t = U^t Z_i psi_0
function ZZ_i(phi_t, psi_t, i::Int, sites::Vector{<:Index})

    ZZ(phi_t, psi_t, i, sites)
end


#----- special case: psi_0 a Z-basis product state (neel_state, random_bitstring)
#
# Z_i psi_0 = z_i psi_0 with z_i=±1, so phi_t = z_i psi_t and
# <Z_i(0) Z_j(t)> = z_i <psi_t|Z_j|psi_t>: only one evolution needed, result real

# <psi|Z_j|psi>, dense ITensor
function z_expect(psi, j::Int, sites::Vector{<:Index})

    real((dag(psi)*apply_Z(psi, j, sites))[])
end

# <psi|Z_j|psi>, MPS
function z_expect(psi::MPS, j::Int, sites::Vector{<:Index})

    real(expect(psi, "Z"; sites=j))
end


function ZZ_zbasis(psi_0, psi_t, i::Int, j::Int, sites::Vector{<:Index})

    z_i=z_expect(psi_0, i, sites)
    abs(abs(z_i)-1.0)<1e-10 || error("ZZ_zbasis: psi_0 is not a Z eigenstate on site $i (<Z_$i> = $z_i), use ZZ with phi_t = U^t Z_i psi_0")

    round(z_i)*z_expect(psi_t, j, sites)
end

# same site: <Z_i(0) Z_i(t)>
function ZZ_i_zbasis(psi_0, psi_t, i::Int, sites::Vector{<:Index})

    ZZ_zbasis(psi_0, psi_t, i, i, sites)
end

#-------------------------------------
# pre compiler
#-------------------------------------

function halfchainee_compile()

L=6;
epsilon=0.9;
thetamean=0.0;


# Initialise the product state and the gates

sites = siteinds("S=1/2", L);
dummysites = siteinds("S=1/2", L);
psi = productMPS(sites, fill("↑", L));

Psi=psi[1]*psi[2]
for i in 3:L
    Psi=Psi*psi[i]
end

brick=brickwall_tensor(L,thetamean,epsilon);
brick=itensorise(brick,sites,dummysites);

e_tev=Array{Float64}(undef, 0)

for t in 1:100
    Psi=brickwall_tev(Psi, brick, sites, dummysites);
    push!(e_tev,halfchainee(Psi,sites))
end

end
   