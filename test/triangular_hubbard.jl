using FiniteLattices

# Independent canonical occupation-basis reference for the triangular tests.
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

@testset "Triangular geometry and physical spin operators" begin
    for (L,W,bx,by,nbonds) in ((2,2,:OBC,:OBC,5),(2,2,:OBC,:PBC,6),
            (2,2,:PBC,:PBC,6),(3,3,:PBC,:PBC,27),(12,4,:OBC,:PBC,136))
        c=tiny_hubbard_config(;L,W,boundary_x=bx,boundary_y=by,symmetry=:u1su2)
        model=hubbard_model(c;lattice_kind=:triangular)
        @test length(model.bonds)==nbonds
        @test allunique(model.bonds) && all(i<j for (i,j) in model.bonds)
        @test model.sublattice===nothing
        if L==12
            bulk=Set(i for i in 1:L*W if model.lattice[i][1] in 5:8)
            @test length(bulk)==16
            @test count(ij->ij[1] in bulk && ij[2] in bulk,model.bonds)==40
        end
    end
    @test_throws ArgumentError hubbard_model(tiny_hubbard_config();lattice_kind=:triangular)
    @test_throws ArgumentError hubbard_spin_correlation_mpo(2,[(1,2),(2,1)])
    p=U1SU2Fermion.n-2U1SU2Fermion.nd
    @test norm(p*p-p)<1e-14
    @test sort(eigvals(Hermitian(convert(Array,p)))) ≈ [0.,0.,1.,1.]
    sl,sr=U1SU2Fermion.SS
    @tensor ss[a b;c d] := sl[a;c v]*sr[v b;d]
    @test eigvals(Hermitian(reshape(convert(Array,ss),16,16))) ≈
        vcat(-.75,zeros(12),fill(.25,3)) atol=1e-13
end

@testset "Doped triangular Hubbard canonical thermal/Born ED" begin
    # Physical square of primitive coordinates, snake indices (1,2,4,3).
    # Independent list has two triangles and five undirected NN bonds.
    bonds=[(1,2),(1,4),(2,4),(2,3),(3,4)]
    for Ne in (3,5)
        # Full four-site purification bond dimension.
        c=tiny_hubbard_config(;L=2,W=2,Ne,U=12.,symmetry=:u1su2,ensemble=:canonical,
            settn_bond_dimension=256,purification_bond_dimension=256,sample_bond_dimension=64)
        model=hubbard_model(c;lattice_kind=:triangular)
        @test Set(model.bonds)==Set(bonds)
        spin=hubbard_spin_correlation_mpo(4,model.bonds)
        reference=TriangularReference.operators(4,Ne,bonds)
        thermal=u1su2_thermal_trajectory(model.H,model.N,c;before_step=i->Random.seed!(4100+i))
        rho=thermal.thermal_mpo
        X=[collect(x) for x in Iterators.product(fill(1:3,4)...) if sum(x.-1)==Ne]
        born=[hubbard_born(rho,c;configuration=x,rng=MersenneTwister(1)) for x in X]
        particles=[hubbard_sample(b,model,c,i) for (i,b) in enumerate(born)]
        @test sum(exp(p.log_p0) for p in particles) ≈ 1 atol=2e-10
        logZ0=last(thermal.log_normalization)+log(born[1].target_sector_weight)+last(thermal.alpha)*Ne
        function step!(state,H,delta,D,index,generation,lineage)
            Random.seed!(4200+lineage)
            env=Environment(state',H,state;disk=false)
            TDVPSweep1!(env,delta;K=32,tol=1e-12,CBEAlg=cbe_algorithm(D,c),
                trunc=tensor_truncation(D,1e-15),GCstep=false,GCsweep=false,verbose=0)
        end
        for (index,beta) in enumerate(sample_betas(c))
            index>1 && foreach(p->advance_hubbard_sample!(p,model,c,index;step!),particles)
            weights=[p.log_p0+p.a for p in particles] # Born-probability-weighted enumeration
            q=logweighted_mean([hubbard_expectation(p.state,spin;disk=false) for p in particles],weights)
            values=[hubbard_local_observables(p.state;disk=false) for p in particles]
            n=logweighted_mean([sum(v.number) for v in values],weights)
            d=logweighted_mean([sum(v.double) for v in values],weights)
            single=(n-2d)/4
            ed=TriangularReference.thermal(reference,beta)
            @test all(p.N==Ne for p in particles)
            @test n ≈ Ne atol=2e-9
            @test logweighted_mean([p.energy for p in particles],weights) ≈ ed.energy atol=8e-7
            @test q ≈ ed.spin_sum atol=8e-8
            @test d ≈ ed.double atol=8e-8
            @test single ≈ ed.single atol=8e-8
            @test 4q/(9*4*single^2) ≈ ed.css_tot atol=8e-8
            @test ed.spin_sum ≈ 3ed.zz atol=1e-12
            @test logZ0+logsumexp(weights) ≈ ed.logZ atol=8e-7
        end
        hot=TriangularReference.thermal(reference,0.)
        @test hot.css_tot ≈ 0 atol=1e-14
        @test hot.single ≈ Ne*(8-Ne)/(4*7) atol=1e-14
    end
end
