
#-------------------------------------
# half chain EntanglementEntropy for itensor
#-------------------------------------

function halfchainee(Psi, sites)

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
   