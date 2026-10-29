
######################################################
# ITensor takes over
######################################################


################## applying the polynomial on a vector ####################

function apply_polynomial(psi, brick, sites, dummysites, m, phi; method=:ITensor)
    psi1=psi
    for i in 1:m
        psi1=psi+exp(-im*phi)*brickwall_tev(psi1, brick, sites, dummysites)
    end
    return psi1
end




# Function to perform Arnoldi iteration
function arnoldi(brick, sites, dummysites, phi, k::Int, m::Int; method=:ITensor)   
    
    v = random_vector(sites)
    # Random initial vector
    H = zeros(ComplexF64, k, k)  # Hessenberg matrix
    V = []  # Orthonormal basis
    
    push!(V, v / norm(v))  # Normalize the initial vector
    
    for j in 1:k-1
        w = apply_polynomial(V[j], brick, sites, dummysites, m, phi; method=:ITensor)  # Apply matrix A to the j-th basis vector
        
        for i in 1:j
            H[i, j] = dot(V[i], w)  # Compute the j-th column of H
        end

        for i in 1:j     
            w -= H[i, j] * V[i]  # Orthogonalize w against previous basis vectors
        end
        
        H[j + 1, j] = norm(w)  # Compute the (j+1)-th entry of the j-th column of H
        
        if H[j + 1, j] != 0
           push!(V,  w / H[j + 1, j])  # Normalize the (j+1)-th basis vector
        end
    end
    
    return H, V
end

function multiply_MpS_mat(V, B)
    N = length(V)
    # Ensure B is of size NxN
    @assert size(B) == (N, N) "Matrix B should be NxN"
    
    # Initialize the result array to store the linear combinations of ITensors
    result = Vector{ITensor}(undef, N)
    
    for i in 1:N
        result[i] = zero(V[1])  # Initialize with a zero ITensor
        for j in 1:N
            result[i] += B[j, i] * V[j]  # Linear combination of ITensors
        end
    end
    
    return result
end


# Function to perform Arnoldi diagonalization
function arnoldi_diagonalize(brick, k::Int, m::Int, phi ,sites, dummysites; method=:ITensor)

    # doing the arnoldi
    H, V = arnoldi(brick, sites, dummysites, phi, k::Int, m::Int ; method=:ITensor)
    # Perform Arnoldi iteration
    
    # Diagonalize H
    evals, evecs = eigen(H[1:k, 1:k])
    
    # Approximate eigenvectors of A
    approx_eigenvectors =  multiply_MpS_mat(V, evecs)
    
    return evals, approx_eigenvectors
end
function extract_eigenvalues(brick, V, sites, dummysites)
    k = length(V)
    eigenvalue_vector_pairs = []

    # Compute eigenvalues and store with corresponding vectors
    for i in 1:k
        v = V[i]
        Av = brickwall_tev(v, brick, sites, dummysites)
        eigenvalue = dot(v, Av) / dot(v, v)  # Rayleigh quotient
        push!(eigenvalue_vector_pairs, (eigenvalue, v))
    end

    # Sort pairs by the absolute value of the eigenvalues, in descending order
    eigenvalue_vector_pairs = sort(eigenvalue_vector_pairs, by = x -> abs(x[1]), rev = true)

    # Keep only the largest half of the eigenvalue-vector pairs
    eigenvalue_vector_pairs = eigenvalue_vector_pairs[1:div(k, 2)]

    # Extract sorted eigenvalues and associated vectors
    eigenvalues = [pair[1] for pair in eigenvalue_vector_pairs]
    sorted_vectors = [pair[2] for pair in eigenvalue_vector_pairs]

    return eigenvalues, sorted_vectors
end

########### Final Wrapper function ####################

function polfed(L,theta,epsilon, k::Int, m::Int, phi;method=:ITensor)
    sites = siteinds("S=1/2", L);
    dummysites = siteinds("S=1/2", L);

		#psi = productMPS(sites,"↑"  )
		#psi = Neel_state(sites  )
		psi = randomMPS(sites  )
		#psi = random_BitString(sites  )

		# Defining Psi on ITensor sites
Psi=psi[1]*psi[2]
for i in 3:L
    Psi=Psi*psi[i]
end

########## Defining Brickwall ####################

brick=brickwall_tensor(L, theta, epsilon)
brick=itensorise(brick, sites, dummysites);

phi=1.0
evals, evecs = arnoldi_diagonalize(brick,  k, m, phi ,sites, dummysites; method=:ITensor);
egivals, eigvecs=extract_eigenvalues(brick, evecs, sites, dummysites);

return egivals, eigvecs
end
