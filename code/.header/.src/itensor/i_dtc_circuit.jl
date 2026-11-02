

function circuit_dtc(sites::Vector{<:Index}, thetamean, epsilon)
    
  #######################
  # Constructing the Brickwall with ITensor
  ########################

    # Number of sites
    N=length(sites)
    # Loading an array of 'site' indices
    s=sites
    # Make gates (1,2),(2,3),(3,4),...

    
  
  hz=rand(N)*2pi;
  Jzz=rand(N)*pi;
  thetadev=pi/50;
  theta=thetamean+randn(1)[]*thetadev;                # Interaction + tuning noise, one draw per circuit as in brickwall_dtc

  # Imperfection in Z tuning, drawn bond by bond (1..N-1) in the same order as brickwall_dtc
  delhs=[randn(4)*pi/50 for j in 1:N-1];

  g=pi*(1-epsilon);
  gatesEven = ITensor[]
  gatesOdd = ITensor[]


  for j in 1:2:(N - 1)
    s1 = s[j]
    s2 = s[j + 1]

    ##### two body gates (Pauli operators, as in brickwall_dtc)
    hj =
      Jzz[j]*op("Z", s1) * op("Z", s2) +
      theta / 2 * op("X", s1) * op("X", s2) +
      theta / 2 * op("Y", s1) * op("Y", s2)

    #### error in two body gates: delh[1],delh[2] before, delh[3],delh[4] after

    delh=delhs[j];

    noise1=exp(-im*delh[3]*op("Sz", s1)) ;
    noise2=exp(-im*delh[4]*op("Sz", s2)) ;


    noise3=exp(-im*delh[1]*op("Sz", s1)) ;
    noise4=exp(-im*delh[2]*op("Sz", s2)) ;

    Gj = exp(-im * hj)
    Gj=apply(noise1,Gj)
    Gj=apply(noise2,Gj)

    Gj=apply(Gj,noise3)
    Gj=apply(Gj,noise4)

    #applying background disorder

    hfildj1=exp(-im*hz[j]*op("Sz", s1)) ;
    hfildj2=exp(-im*hz[j+1]*op("Sz", s2)) ;

    Gj=apply(Gj,hfildj1)
    Gj=apply(Gj,hfildj2)
    
    push!(gatesEven, Gj)
  end

   

  for j in 2:2:(N - 2)
    s1 = s[j]
    s2 = s[j + 1]

    ##### two body gates (Pauli operators, as in brickwall_dtc)
    hj =
      Jzz[j]*op("Z", s1) * op("Z", s2) +
      theta / 2 * op("X", s1) * op("X", s2) +
      theta / 2 * op("Y", s1) * op("Y", s2)

    #### error in two body gates: delh[1],delh[2] before, delh[3],delh[4] after

    delh=delhs[j];

    noise1=exp(-im*delh[3]*op("Sz", s1)) ;
    noise2=exp(-im*delh[4]*op("Sz", s2)) ;


    noise3=exp(-im*delh[1]*op("Sz", s1)) ;
    noise4=exp(-im*delh[2]*op("Sz", s2)) ;

    Gj = exp(-im * hj)
    Gj=apply(noise1,Gj)
    Gj=apply(noise2,Gj)

    Gj=apply(Gj,noise3)
    Gj=apply(Gj,noise4)
    
    #### adding Kicks

  
    xj1=exp(-im*g*op("Sx", s1))
    xj2=exp(-im*g*op("Sx", s2))

    Gj=apply(xj1,Gj)
    Gj=apply(xj2,Gj)


    push!(gatesOdd, Gj)
  end


  #### Boundary Kicks

  s1 = s[1]
  sN = s[N]

  xj1=exp(-im*g*op("Sx", s1))
  xjN=exp(-im*g*op("Sx", sN))


  gatesEven[1]=apply(xj1,gatesEven[1])
  gatesEven[end]=apply(xjN,gatesEven[end])




      
  return gatesEven, gatesOdd , Jzz, hz
end

