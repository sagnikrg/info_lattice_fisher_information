#######################################
# The Krylov part of the power basis on a GPU: same algorithm as arnoldi_power (functions/power_basis.jl).
# Not included by polfed.jl, so nothing else depends on CUDA: after functions/overlap.jl and power_basis.jl,
#   using CUDA
#   include(".header/.src/polfed/functions/power_basis_gpu.jl")
#######################################

using CUDA
using CUDA.CUFFT: plan_fft!, plan_ifft!


#####################################
# S x with the circulant embedding of S, on the GPU
#####################################

# S_ij = c_{j-i}, its first column conj(c), its first row c; embedded into a circulant matrix of size 2n

struct ToeplitzGPU{TC,TP,TI}
    n::Int
    fcol::TC                     # the transform of the first column of the circulant matrix
    buf::TC                      # work space of length 2n
    fwd::TP
    inv::TI
end

function ToeplitzGPU(columneinz::AbstractVector, n::Int)
    col=zeros(ComplexF64, 2n)
    col[1:n]=conj.(columneinz[1:n])
    col[n+2:2n]=reverse(columneinz[2:n])
    buf=CuArray(col)
    fwd=plan_fft!(buf); inv=plan_ifft!(buf)
    fwd*buf
    return ToeplitzGPU(n, buf, CUDA.zeros(ComplexF64, 2n), fwd, inv)
end

function toeplitz_mul!(y::AbstractVector, F::ToeplitzGPU, x::AbstractVector)
    fill!(F.buf, 0)
    view(F.buf, 1:length(x)).=x
    F.fwd*F.buf
    F.buf.*=F.fcol
    F.inv*F.buf
    y.=view(F.buf, 1:F.n)
    return y
end


#####################################
# Arnoldi: V, SV on the GPU, H on the CPU; g_k is a running sum, serial, so it stays on the CPU
#####################################

function arnoldi_power_gpu(columneinz::AbstractVector, k::Int, ncv::Int, phi::Real)

    n=ncv*k + 1
    F=ToeplitzGPU(columneinz, n)
    V=CUDA.zeros(ComplexF64, n, ncv); SV=CUDA.zeros(ComplexF64, n, ncv)
    H=zeros(ComplexF64, ncv+1, ncv)
    w=CUDA.zeros(ComplexF64, n); Sw=CUDA.zeros(ComplexF64, n)
    v=zeros(ComplexF64, n)

    v[1]=1
    copyto!(view(V, :, 1), v)
    toeplitz_mul!(view(SV, :, 1), F, view(V, :, 1))

    for j in 1:ncv
        len=(j-1)*k + 1; lw=min(len + k, n)
        fill!(w, 0)
        a=Array(view(V, 1:len, j))
        copyto!(view(w, 1:lw), apply_polynomial(a, k, phi)[1:lw])

        vw=view(w, 1:lw); Vj=view(V, 1:lw, 1:j); SVj=view(SV, 1:lw, 1:j)
        for pass in 1:2
            h=SVj'*vw
            H[1:j, j]+=Array(h)
            vw.-=Vj*h
        end

        toeplitz_mul!(Sw, F, w)
        H[j+1, j]=sqrt(max(real(dot(vw, view(Sw, 1:lw))), 0.0))
        if j < ncv
            if abs(H[j+1, j]) < 1e-14                    # the end of the Krylov space (Chapter 7)
                return H[1:j, 1:j], Array(view(V, :, 1:j))
            end
            view(V, :, j+1).=w ./ H[j+1, j]
            view(SV, :, j+1).=Sw ./ H[j+1, j]
        end
    end

    return H[1:ncv, 1:ncv], Array(V)
end
