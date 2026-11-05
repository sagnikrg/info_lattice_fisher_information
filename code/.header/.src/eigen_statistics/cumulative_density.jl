## functions for cumulative density of level spacing statistics ##

### CDF(r) = probability of a ratio less than r, without any binning: a step of 1/N at each of the N ratios.
### r is always the flat list of ratios, e.g. LevelSpacingRatio(E) of the nested output of the symmetries functions.
### For several realisations pool the ratios first, vcat(r1, r2, ...)

include("../random_matrices/random_matrices.jl")


################################
# CDF of the data
################################


# The step function itself, returns the sorted ratios and the CDF at each of them, for plotting

function CumulativeDensity(r)
    x=sort(r)
    return x, (1:length(x))/length(x)
end


# The step function evaluated on the points x

function CumulativeDensity(r, x)
    rsorted=sort(r)
    return [searchsortedlast(rsorted, xi) for xi in x]/length(rsorted)
end



################################
# Reference curves
################################


CumulativeDensityPoisson(x) = 2 .- 2 ./ (x .+ 1)


# Ratios of a random matrix ensemble: :CUE, :COE or :GOE
# for :GOE only the central 80% of the spectrum is used

function RandomMatrixRatios(ensemble::Symbol; n=200, nsamples=500)

    r=Float64[]

    for _ in 1:nsamples

        if ensemble==:CUE
            append!(r, LevelSpacingRatio(eigvals(CUE(n))))
        elseif ensemble==:COE
            append!(r, LevelSpacingRatio(eigvals(COE(n))))
        elseif ensemble==:GOE
            E=eigvals(Symmetric(GOE(n)))
            cut=round(Int, n/10)
            E=E[cut+1:n-cut]
            append!(r, LevelSpacingRatio(cis.(pi*E/maximum(abs.(E))*0.99)))   # LevelSpacingRatio takes phases
        else
            error("Unknown ensemble $ensemble")
        end

    end

    return r
end



################################
# Difference to a reference
################################


# Delta CDF on the points x, data minus reference
# reference is a function of x, like CumulativeDensityPoisson, or a list of ratios, like RandomMatrixRatios(:COE)

CumulativeDensityDifference(r, reference::Function, x) = CumulativeDensity(r, x) .- reference(x)

CumulativeDensityDifference(r, reference::AbstractVector, x) = CumulativeDensity(r, x) .- CumulativeDensity(reference, x)
