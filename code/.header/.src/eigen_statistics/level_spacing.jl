### functions for level spacing and level spacing ratio 
###
### including one for directly generating the histogram plot


### functions for building eigenstatics

using StatsBase


################################
# Level Spacing Statistics
################################ 


# Exact degeneracies: levels closer than tol are one level, each distinct level is counted once.
# Otherwise two consecutive zero spacings make the ratio 0/0. The angles -pi and pi are also the same level.
# tol is far below the mean spacing 2pi/2^L for all L we can diagonalise, and above the roundoff of eig.

function DistinctLevels(EigA; tol=1e-10)

    N=sort(angle.(EigA),rev=true);                  # descending order
    length(N)<2 && return N

    N=N[[true; (N[1:end-1]-N[2:end]).>tol]];        # the first one of every group
    length(N)>1 && N[1]-N[end]>2*pi-tol && pop!(N);

    return N
end


function LevelSpacing(EigA; tol=1e-10)

    N=DistinctLevels(EigA; tol=tol);   # N is the array of distinct eigenvalues in angle form in descending order
    length(N)<2 && return Float64[]

    EigA1=copy(N);                  # EigA1 is the array of eigenvalues in angle form in descending order
   
    ls = deleteat!(EigA1,1);        # ls is the array of eigenvalues in angle form in descending order without the first element
    la = deleteat!(N,length(N));    # la is the array of eigenvalues in angle form in descending order without the last element            
    
    m=copy(la-ls)/mean(la-ls);      # m is the array of level spacing statistics

    m

end    


function LevelSpacingRatio(EigA; tol=1e-10)

    N=DistinctLevels(EigA; tol=tol);
    length(N)<3 && return Float64[]


    EigA1=copy(N);
    ls = deleteat!(EigA1,1);
    la = deleteat!(N,length(N));
    m=copy(la-ls)/mean(la-ls)

    #histogram(m)

    #Computing Level Spacing Ratio

    n=fill(0.1,length(m)-1);

        for i in 1:length(n)
            #n[i]=m[i+1]/m[i];
            n[i]= minimum([m[i], m[i+1]])/maximum([m[i], m[i+1]]);
        end

    #histogram(n)
    #mean(n)
    n
end    


function MeanLevelSpacingRatio(EigA; tol=1e-10)
    return mean(LevelSpacingRatio(EigA; tol=tol))
end



################################
# Sector Resolved
################################

# EigA is the nested output of the symmetries functions, [[eigenvalues of sector 1], [eigenvalues of sector 2], ...]
# The statistic is computed inside every sector, and the results are pooled.
# MeanLevelSpacingRatio then is the mean over all the ratios, i.e. sectors are weighted by their size.
# Sectors with less than 3 levels have no ratio and are skipped.
# Per sector values: LevelSpacingRatio.(EigA)
# Exact degeneracies are treated inside every sector, as above.


function LevelSpacing(EigA::AbstractVector{<:AbstractVector}; tol=1e-10)
    return reduce(vcat, [LevelSpacing(E; tol=tol) for E in EigA if length(E)>=3]; init=Float64[])
end


function LevelSpacingRatio(EigA::AbstractVector{<:AbstractVector}; tol=1e-10)
    return reduce(vcat, [LevelSpacingRatio(E; tol=tol) for E in EigA if length(E)>=3]; init=Float64[])
end