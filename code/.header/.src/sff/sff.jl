# functions related to spectral form factors

function LevelSpacingExtended(EigA)

    En=angle.(EigA);
    Em=copy(En);
        
    Eigdiff=fill(0.0*im, 2^L, 2^L);
        
    for j in 1:2^L
        for k in 1:2^L
            Eigdiff[j,k]=En[j]-Em[k];
        end
    end
    
    m    

end



################################
# SFF 
################################

function sff(En,t)
    Em=copy(En);
		  
    Eigdiff=fill(0.0*im, length(En), length(En));

    for j in 1:length(En)
        for k in 1:length(En)
            Eigdiff[j,k]=(En[j]/Em[k])^t;
        end
    end

    return real(sum(Eigdiff))
end
