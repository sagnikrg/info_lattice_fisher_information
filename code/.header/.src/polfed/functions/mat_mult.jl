
include("overlap.jl")


################################
# Overlap times Abstract vector
################################

function mat_mult(A::AbstractMatrix, B::AbstractVector)
  
    dimB=length(B)

    # restrict the matrix A to the size of B

    A = A[1:dimB, 1:dimB]

    # Perform matrix-vector multiplication
    C = A*B
    # Check if the result is a vector

    return C
end


########################################
# Restrited inner product
########################################

function inner_prod(vi::AbstractVector ,vj::AbstractVector)

    m= min(length(vi), length(vj))

    # Inner product of two vectors

    return sum(conj(vi[1:m]) .* vj[1:m])

end
 