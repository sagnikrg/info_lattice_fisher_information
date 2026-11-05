### The dtc circuit without the noise, for looking at exact symmetries
###
###     U = K Wodd_even Wodd D ,  K = exp(-i g/2 sum_j X_j) , g = pi(1-epsilon) , D = exp(-i/2 sum_j h_j Z_j)
###     gates exp(-i J_j Z_j Z_{j+1} - i theta/2 (X_j X_{j+1} + Y_j Y_{j+1})) , odd bonds act first
###
### Same layers and same order as circuit_dtc, minus delh and the theta noise.
### Works for any L, every site is kicked exactly once.

using LinearAlgebra


function circuit_dtc_clean(L::Int, theta, epsilon, h, J)

    Paulis=[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]

    g=pi*(1-epsilon)
    site(A,j)=kron(Matrix{ComplexF64}(I,2^(j-1),2^(j-1)), A, Matrix{ComplexF64}(I,2^(L-j),2^(L-j)))
    bond(A,j)=kron(Matrix{ComplexF64}(I,2^(j-1),2^(j-1)), A, Matrix{ComplexF64}(I,2^(L-j-1),2^(L-j-1)))

    Xs,Ys,Zs = Paulis[2],Paulis[3],Paulis[4]

    gate(j)=exp(Matrix{ComplexF64}(-im*J[j]*kron(Zs,Zs)-im*theta/2*(kron(Xs,Xs)+kron(Ys,Ys))))

    U=Matrix{ComplexF64}(I,2^L,2^L)

    for j in 1:L                                                 # the fields
        U=site(exp(Matrix{ComplexF64}(-im*h[j]/2*Zs)),j)*U
    end
    for j in 1:2:L-1                                             # odd bonds
        U=bond(gate(j),j)*U
    end
    for j in 2:2:L-1                                             # even bonds
        U=bond(gate(j),j)*U
    end
    for j in 1:L                                                 # the kick
        U=site(exp(Matrix{ComplexF64}(-im*g/2*Xs)),j)*U
    end

    return U
end


# implicitly drawing the disorder, same intervals as circuit_dtc

circuit_dtc_clean(L::Int, theta, epsilon) = circuit_dtc_clean(L, theta, epsilon, rand(L)*2*pi, rand(L-1)*pi)
