#####################################################################
# Gullans-Huse reference qubit probe on the pi mode (K.A.R.E.N., Ch5)
# R starts in a Bell pair with site x0 (edge x0=1 or bulk x0=L/2), the
# other sites in a random Z product state, and the Floquet unitary U
# (any dense U, site 1 most significant) acts on the chain alone.
# Late time signals over t1..t2:
#   I_R     = mean_t I(R:x0)                       (entropic, blind to 0 vs pi)
#   Api_ZZ  = mean_t (-1)^t <Z_R Z_x0>,  A0_ZZ = mean_t <Z_R Z_x0>
#   (same for XX and YY)
# <Z_R Z_x0>(t) = <Z_x0(t) Z_x0>_config, so the ZZ pi signal at x0=1 is the
# edge pi mode seen by a joint measurement of R and the edge qubit.
#####################################################################
using LinearAlgebra
include(joinpath(@__DIR__, "..", "gullans_huse", "gh_probe.jl"))

function gh_edge_run(U::AbstractMatrix, L::Int, x0::Int, config; t1 = 100, t2 = 200)
    Ut = transpose(U)
    M = gh_initial_state(L, x0, config)
    IR = 0.0; S = Dict(:ZZ => [0.0, 0.0], :XX => [0.0, 0.0], :YY => [0.0, 0.0]); n = 0
    for t in 0:t2
        if t >= t1
            IR += mutual_info_R(M, [x0], L); n += 1
            for (k, P) in ((:ZZ, PZ), (:XX, PX), (:YY, PY))
                c = R_site_correlator(M, P, x0, L)
                S[k][1] += (-1)^t * c; S[k][2] += c
            end
        end
        t < t2 && (M = M * Ut)
    end
    return (I_R = IR / n, Api_ZZ = S[:ZZ][1] / n, A0_ZZ = S[:ZZ][2] / n,
            Api_XX = S[:XX][1] / n, A0_XX = S[:XX][2] / n, Api_YY = S[:YY][1] / n, A0_YY = S[:YY][2] / n)
end
