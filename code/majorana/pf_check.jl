include("majorana_floquet.jl")

function pf_brute(A)
    n=size(A,1); n==0 && return 1.0
    s=0.0
    for j in 2:n
        idx=setdiff(1:n,[1,j])
        s+= (-1)^(j) * A[1,j]*pf_brute(A[idx,idx])
    end
    s
end
ok=true
for n in (2,4,6,8), _ in 1:200
    A=randn(n,n); A=A-A'
    global ok &= sign(pf_brute(A))==pfaffian_sign(A)
end
println("pfaffian sign ok: ",ok)
