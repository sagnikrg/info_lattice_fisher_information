 using ToeplitzMatrices



######################################
# Toeplitz matrix for the overlap matrix
######################################

# columneinz[j+1] = c_j = <E|U^j|E>, j = 0, 1, 2, ...
# In the power basis |j> = U^j|E> the overlap matrix is
#
#   overlap_m[i, j] = <i|j> = <E|U^{j-i}|E> = c_{j-i},   c_{-n} = conj(c_n)
#
# a Hermitian Toeplitz matrix, fixed by its first row (c_0, c_1, c_2, ...) and its first column
# (c_0, c_1^*, c_2^*, ...). c_0 = <E|E> = 1 has to be real.

function overlap(columneinz)

    first_column=conj.(columneinz)
    first_row=columneinz

    ### format <i|j> = c_{j-i}

    overlap=Toeplitz(first_column, first_row)
    return overlap
end


# the same shifted by one step, overlap_shift[i, j] = <i|U|j> = c_{j+1-i}
# (n x n, needs c_0 ... c_n)

function overlap_shift(columneinz, n::Int)

    c=columneinz
    first_row=c[2:n+1]
    first_column=vcat(c[2], c[1], conj.(c[2:n-1]))

    return Toeplitz(first_column, first_row)
end


######################################
# Factorised Toeplitz matrices, for many products with the same matrix
######################################

# S*x with the Toeplitz type plans the FFT (FFTW) and transforms the embedding of S again for every
# product. factorize(S) does both once (2 FFTs per product instead of 3, and no planning), so it is
# kept for the whole Krylov part. x may be shorter than S (zero padded).

overlap_factorized(columneinz, n::Int)=factorize(overlap(columneinz[1:n]))
overlap_shift_factorized(columneinz, n::Int)=factorize(overlap_shift(columneinz, n))

toeplitz_mul!(y::AbstractVector, F, x::AbstractVector)=mul!(y, F, x, true, false)
toeplitz_mul(F, x::AbstractVector, n::Int)=toeplitz_mul!(zeros(ComplexF64, n), F, x)
