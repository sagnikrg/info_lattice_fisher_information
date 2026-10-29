

#----------------------------------------------------------
# homebrew wrappers for tev 
#----------------------------------------------------------

# This file contains the homebrew wrappers for the ITensor library in Julia.

#------------------------------------
# Wrapper for initial random states
#------------------------------------

function random_vector(sites)

    L=length(sites)
		#psi = productMPS(sites,"↑"  )
		#psi = neel_state(sites)
		psi = randomMPS(sites  )
		#psi = random_bitstring(sites  )

		# Defining Psi on ITensor sites
Psi=psi[1]*psi[2]
for i in 3:L
    Psi=Psi*psi[i]
end

return Psi
end


#-------------------------------------
# Néel state
#-------------------------------------

function neel_state(sites)
    # Create a site set for N spin-half sites
   

    # Initialize the MPS as a Néel state
    state = productMPS(sites, n -> isodd(n) ? "↑" : "↓")

    return state
end


#------------------------------------
# Random Bit String
#------------------------------------

function random_bitstring(sites)
    # Create a site set for N spin-half sites

    N=length(sites)
  # Create an array to store the random states
random_states = Vector{String}(undef, N)

# Fill the array with random "↑" or "↓" states
for i in 1:N
    random_states[i] = rand(Bool) ? "↑" : "↓"
end

# Create the MPS representing the product state with the random bit string
psi = productMPS(sites, random_states)

    return psi
end

