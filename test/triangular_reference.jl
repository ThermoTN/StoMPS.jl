# Independent fixed-N occupation-basis reference for the triangular tests.
module TriangularReference
using LinearAlgebra

function hop(state, create, annihilate)
    a, b = 1 << (create-1), 1 << (annihilate-1)
    iszero(state & b) && return nothing
    sign = isodd(count_ones(state & (b-1))) ? -1.0 : 1.0
    state = xor(state,b)
    !iszero(state & a) && return nothing
    sign *= isodd(count_ones(state & (a-1))) ? -1.0 : 1.0
    return state | a, sign
end

function operators(nsites, Ne, bonds; t=1.0, U=12.0,
        bulk_sites=collect(1:nsites), spin_bonds=bonds)
    basis = [s for s in 0:(1 << (2nsites))-1 if count_ones(s)==Ne]
    lookup = Dict(s=>i for (i,s) in enumerate(basis))
    H = zeros(length(basis),length(basis)); Q = similar(H); fill!(Q,0)
    zz, number, double = zeros(length(basis)), zeros(length(basis)), zeros(length(basis))
    for (column,s) in enumerate(basis)
        up = [Int(!iszero(s & (1 << (2i-2)))) for i in 1:nsites]
        down = [Int(!iszero(s & (1 << (2i-1)))) for i in 1:nsites]
        z = (up .- down)./2
        H[column,column] = U*sum(up.*down)
        number[column] = sum(up[i]+down[i] for i in bulk_sites)
        double[column] = sum(up[i]*down[i] for i in bulk_sites)
        for (i,j) in bonds, spin in 0:1, (a,b) in ((i,j),(j,i))
            moved = hop(s,2a-1+spin,2b-1+spin)
            moved===nothing && continue
            H[lookup[moved[1]],column] -= t*moved[2]
        end
        for (i,j) in spin_bonds
            zz[column] += z[i]*z[j]
            # S_i^+ S_j^- and S_i^- S_j^+, applied right to left.
            for (a,b,c,d) in ((2i-1,2i,2j,2j-1),(2i,2i-1,2j-1,2j))
                first = hop(s,c,d)
                first===nothing && continue
                second = hop(first[1],a,b)
                second===nothing && continue
                Q[lookup[second[1]],column] += first[2]*second[2]/2
            end
        end
        Q[column,column] += zz[column]
    end
    return (;H,Q,zz=Diagonal(zz),number=Diagonal(number),double=Diagonal(double),
        single=Diagonal((number.-2double)./length(bulk_sites)),Nbulk=length(bulk_sites))
end

function thermal(ops,beta)
    rho = exp(-beta*Hermitian(ops.H))
    Z = tr(rho)
    mean(op) = real(tr(rho*op)/Z)
    q,p = mean(ops.Q),mean(ops.single)
    return (;energy=mean(ops.H),spin_sum=q,zz=mean(ops.zz),number=mean(ops.number),
        double=mean(ops.double),single=p,logZ=log(Z),css_tot=4q/(9ops.Nbulk*p^2))
end
end
