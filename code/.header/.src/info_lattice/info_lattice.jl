###################################################
# functions to compute information lattice
####################################################



function entanglement_entropy(eigenstate::AbstractVector)
    # Compute the density matrix of the eigenstate
    density_matrix = eigenstate * eigenstate'

    # Compute the eigenvalues of the density matrix
    eigenvalues = eigvals(density_matrix)

    # Compute the entanglement entropy
    entropy = -sum(eigenvalues .* log.(eigenvalues))

    return entropy
end

function replace_small_values(arr::Vector{Float64})
    new_arr = similar(arr)  # Create a new array of the same size and type
    for i in eachindex(arr)
        if arr[i] < 1e-15
            new_arr[i] = 1e-15
        else
            new_arr[i] = arr[i]
        end
    end
    return new_arr
end


function von_neumann_entropy(eigenvalues::AbstractVector; base::Real = 2.0)  # add cut off
    # Filter out zero eigenvalues to avoid log(0) errors
    valid_eigenvalues = replace_small_values(real.(eigenvalues))
   
    # Compute entropy
    entropy = -sum(λ -> λ * log(λ) / log(base), valid_eigenvalues)
    return entropy
end




function site_entropy(Psi, i, l, L)
    
psi=reshape(Psi,(2^(i-1),2^(l),2^(L-i-l+1)))
psi_permuted = permutedims(psi, (1, 3, 2))
psi=reshape(psi_permuted,(2^(L-l),2^l))

rho=psi'*psi
eigvals, eigvecs=eigen(rho)
#von_neumann_entropy(eigvals)

return von_neumann_entropy(eigvals)
end





function site_bit(Psi, i, l, L)
    
    psi=reshape(Psi,(2^(i-1),2^(l),2^(L-i-l+1)))
    psi_permuted = permutedims(psi, (1, 3, 2))
    psi=reshape(psi_permuted,(2^(L-l),2^l))
    
    rho=psi'*psi
    eigvals, eigvecs=eigen(rho)
    #von_neumann_entropy(eigvals)
    
    return l-von_neumann_entropy(eigvals)
end




function info_lattice_raw(Psi, L)
    # Predefine the type as a vector of Float64 vectors
    info_lattice_raw_t = Vector{Vector{Float64}}(undef, L)

    for l in 1:L
        bit_info = Float64[]  # Store bit information for sub-lattices of size `l`
        for i in 1:(L - l + 1)
            push!(bit_info, site_bit(Psi, i, l, L))  # Compute and store bit information
        end
        info_lattice_raw_t[l] = bit_info  # Assign explicitly typed row
    end

    return info_lattice_raw_t
end




function info_lattice_bug(Psi, L)
    
    info_lattice_raw_t=info_lattice_raw(Psi,L)
    info_lattice_t=info_lattice_raw_t
    info_lattice_t[1].=info_lattice_raw_t[1]

    bit_info=Float64[]
    
    for i in 1:(L-1)
        push!(bit_info, info_lattice_raw_t[2][i]-info_lattice_raw_t[1][i]-info_lattice_raw_t[1][i+1])
    end

    info_lattice_t[2].=bit_info

    for l in 3:L
        bit_info=Float64[]
        for i in 1:(L-l+1)
            push!(bit_info, info_lattice_raw_t[l][i]-info_lattice_raw_t[l-1][i]-info_lattice_raw_t[l-1][i+1]+info_lattice_raw_t[l-2][i+1])
        end
        info_lattice_t[l].=bit_info
    end

    return info_lattice_t
end



function info_lattice(Psi)
    
    L=size(Psi,1) |> log2 |> Int

    info_lattice_raw_t=info_lattice_raw(Psi,L)
    info_lattice_t= Vector{Vector{Float64}}(undef, L)
    info_lattice_t[1]=info_lattice_raw_t[1]

    bit_info=Float64[]
    
    for i in 1:(L-1)
        push!(bit_info, info_lattice_raw_t[2][i]-info_lattice_raw_t[1][i]-info_lattice_raw_t[1][i+1])
    end

    info_lattice_t[2]=bit_info

    for l in 3:L
        bit_info=Float64[]
        for i in 1:(L-l+1)
           push!(bit_info, info_lattice_raw_t[l][i]-info_lattice_raw_t[l-1][i]-info_lattice_raw_t[l-1][i+1]+info_lattice_raw_t[l-2][i+1])
        end
        info_lattice_t[l]=bit_info
    end

    return info_lattice_t
end



function info_lattice(U::Matrix{ComplexF64})

 

    eigvals,eigvecs=eigen(U);

    Psi=eigvecs[:,1]
    infolattice=info_lattice(Psi)

    for i in 2:length(eigvals)
        Psi=eigvecs[:,i]
        infolattice+=info_lattice(Psi)
    end

    infolattice=infolattice./length(eigvals)
end




#### aux

function info_lattice_brickwall(L::Int, theta::Float64, epsilon::Float64)

 
    U=circuit_dtc(L,theta, epsilon)
    eigvals,eigvecs=eigen(U);

    Psi=eigvecs[:,1]
    infolattice=info_lattice(Psi)

    for i in 2:length(eigvals)
        Psi=eigvecs[:,i]
        infolattice+=info_lattice(Psi)
    end

    infolattice=infolattice./length(eigvals)
end




function info_lattice_proj(U::Matrix{ComplexF64}, indlist::Array{Int})

 

    eigvals,eigvecs=eigen(U);

    Psi=eigvecs[:,1]
    Psi=projectatt(Psi, indlist)
    infolattice=info_lattice(Psi)

    for i in 2:length(eigvals)
        Psi=eigvecs[:,i]
        Psi=projectatt(Psi, indlist)
        infolattice+=info_lattice(Psi)
    end

    infolattice=infolattice./length(eigvals)
end


include("draw_info_lattice.jl")


function lattice_to_vec(info_lattice_state)
   
    L=size(info_lattice_state)[1]

    A=info_lattice_state[1][:]

    for i in 2:L 
        A=cat(A,info_lattice_state[i][:], dims=1 )
    end

    return A  
end