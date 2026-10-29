

function Brickwall(thetamean, epsilon, sites; method=:ITensor)
    
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
          
  g=pi*(1-epsilon);
  gatesEven = ITensor[]
  gatesOdd = ITensor[]


  for j in 1:2:(N - 1)
    s1 = s[j]
    s2 = s[j + 1]
    theta=thetamean+randn(1)[]*thetadev;

    ##### two body gates
    hj =
      Jzz[j]*op("Sz", s1) * op("Sz", s2) +
      theta / 2 * op("Sx", s1) * op("Sx", s2) +
      theta / 2 * op("Sy", s1) * op("Sy", s2)

    #### error in two body gates

    delh=randn(4)*pi/50;
    
    noise1=exp(-im*delh[1]*op("Sz", s1)) ;
    noise2=exp(-im*delh[2]*op("Sz", s2)) ;

    
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
    theta=thetamean+randn(1)[]*thetadev;

    ##### two body gates
    hj =
      Jzz[j]*op("Sz", s1) * op("Sz", s2) +
      theta / 2 * op("Sx", s1) * op("Sx", s2) +
      theta / 2 * op("Sy", s1) * op("Sy", s2)

    #### error in two body gates

    delh=randn(4)*pi/50;
    
    noise1=exp(-im*delh[1]*op("Sz", s1)) ;
    noise2=exp(-im*delh[2]*op("Sz", s2)) ;

    
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

