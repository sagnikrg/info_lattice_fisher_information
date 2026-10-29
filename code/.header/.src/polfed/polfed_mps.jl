#######################################
# POLFED with linear algebra functions
#######################################

###
# The entire idea of Polfed is based on efficient matrix-vector multiplication.
# Let's say we want the eigenvalues of the polynomial 1+U+U^2....+U^m
# and with k krylov vectors.
# then all of the rest of the operations can be expressed as linear algebra operations on the overlaps <psi_0|U^j\daggerU^i|psi_0>
# where |psi_0> is the initial vector.
# and i, j runs from 0 to mk

# The following function calculates the overlap <psi_0|U^j\daggerU^i|psi_0>

# Now in principle, the whole matrix can be generated just from the first column. For general purpose we will use this as the input to our algorithm.


## just testing gitpush




#####################################
# Efficient generation of inner products
#####################################  

# Recall: one needs psi_0.... psi_{k-1}
# Recall: for each degree m of the Polynomial there is m+1 term

# wishlist: ideally should be able to generate the inner products in a single sweep (added in main body)

using SparseArrays


#####################################
# Now, since we are working in a non-orthogonal basis, we need to calculate the inner products of any two vectors
# expressed in this basis, in terms of the overlaps <psi_0|U^j\daggerU^i|psi_0>

# The following function calculates the inner product of two vectors in the basis of the overlaps <psi_0|U^j\daggerU^i|psi_0>

# The input is the two indices i and j of the vectors, the coefficients of the vectors in the basis of the overlaps, and the overlaps themselves

# The output is the inner product of the two vectors

# The function is defined as follows:

function inner_prod(vi,vj, overlap)
    n_coeff=length(vi)
    sum=0.0*im
    for k in 1:n_coeff
         for l in 1:n_coeff
             sum=sum+conj(vi[k])*vj[l]*overlap[l,k]
         end
     end
     return sum
 end
 


function inner_prod(i,j, coeff, overlap)
   n_coeff=size(coeff)[2]
   sum=0.0*im
   for k in 1:n_coeff
        for l in 1:n_coeff
            sum=sum+conj(coeff[i,k])*coeff[j,l]*overlap[l,k]
        end
    end
    return sum
end


#############
# And we need a way to apply the polynomial on the vector
###########

function apply_polynomial(coeff,l, m, k, phi )
    n_coeff=size(coeff)[2]
    for i in 0:(m)*(l-2)
        for j in 0:m
         
         coeff[l,(j+i+1)]+=(exp(im*phi)^j)*coeff[l-1,i+1]  # this step basically ammounts to applying the polynomial to the previous vector
        end
    end

    return coeff
end

#############
# Now that we have the eigenvector of the Hessenberg matrix, we need to restore the eigenvector of the original matrix
# The following function restores the eigenvectors of the original matrix from the eigenvectors of the Hessenberg matrix


function restore_eigvecs(V, B)
    N = size(V)[1]
    # Ensure B is of size NxN
    @assert size(B) == (N, N) "Matrix B should be NxN"
    
    # Initialize the result array to store the linear combinations of ITensors
    #result = Vector{ITensor}(undef, N)
    
    result=zeros(ComplexF64, size(V))

    for i in 1:N
        for j in 1:N
            result[i,:] .+= B[j, i] * V[j,:]  # Linear combination of ITensors
        end
    end
    
    return result
end

#############
# Finally we need to compute the eigenvalues of the original Unitary
# Note that it is not necessary to compute the eigenvectors of the original Unitary, one can just read it off from the elements <i|j+1>


function eigenvalues(eigvecscoeffs, overlap_m)
    
    eigenvallength=size(eigvecscoeffs,1)

    eigvalsfull=fill(0.0*im,eigenvallength)

    ncoeffs=size(eigvecscoeffs,2)



    for i in 1:eigenvallength
        sum=0.0*im
        sum_norm=0.0*im
            for j in 1:ncoeffs
                for l in 1:ncoeffs
                sum+=conj(eigvecscoeffs[i,j])*eigvecscoeffs[i,l]*overlap_m[l,j+1]
                sum_norm+=conj(eigvecscoeffs[i,j])*eigvecscoeffs[i,l]*overlap_m[l,j]
            end

        end
    eigvalsfull[i]=sum/sum_norm
    end

    return eigvalsfull
end

function convergence(eigvecscoeffs, eigenvalues, overlap_m)
    n=length(eigenvalues)
    convergence_modulus = zeros(ComplexF64,n-1)
    for i in 1:n-1
        psi = eigvecscoeffs[i, :]
        push!(psi, 0.0)
        lambda = eigenvalues[i]
        psi1=zeros(ComplexF64, length(psi))
        for j in 1:length(psi)-1
            psi1[j+1]=psi[j]
        end
        #inner_prod(psi1, psi1, overlap_m)
        convergence_modulus[i] =  + abs.(lambda)^2 * inner_prod(psi, psi, overlap_m) #- conj.(lambda) * inner_prod(psi, psi1, overlap_m)-lambda*inner_prod(psi1, psi, overlap_m)
    end
    return convergence_modulus
end