using ToeplitzMatrices



######################################
# Sparse matrix for the overlap matrix
######################################

function overlap(columneinz)

    first_column=columneinz
    first_row=conj.(columneinz)
    
    ### format <0|1>
    
    overlap=Toeplitz(first_column,first_row)
    return overlap
end
    