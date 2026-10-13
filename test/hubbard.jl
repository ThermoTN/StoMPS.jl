using LinearAlgebra, Random, FiniteMPS

function tiny_hubbard_config(;kwargs...)
    return merge((;L=2,W=1,boundary_x=:OBC,boundary_y=:OBC,t=1.0,U=8.0,mu=4.0,
        symmetry=:su2su2,ensemble=:half_filling,Ne=2,Ns=4,
        beta_list=[1/128,1/16,1/8],beta0=1/16,
        settn_bond_dimension=32,purification_bond_dimension=32,sample_bond_dimension=32,
        cbe_ratio=.25,cbe_tolerance=1e-8,rsvd=true,krylov_dimension=32,tdvp_tolerance=1e-12,
        settn_maxorder=6,settn_maxiter=8,settn_noise=[(.1,.1),(.1,.01),(.1,.001)],
        settn_tolerance=1e-13,settn_compression=1e-16,settn_cutoff=1e-15,
        purification_cutoff=1e-15,born_cutoff=1e-15,sample_cutoff=1e-15,
        probability_tolerance=1e-10,number_tolerance_per_site=1e-9,
        number_flow_maximum_iterations=24,environment_disk=false,gc_step=false,gc_sweep=false),
        (;kwargs...))
end

@testset "Hubbard geometry and thermal/Born ED" begin
    for (L,W,bx,by,expected) in ((4,1,:OBC,:OBC,3),(4,1,:PBC,:OBC,4),
            (2,2,:OBC,:OBC,4),(2,2,:OBC,:PBC,4),(4,2,:PBC,:PBC,12))
        model=hubbard_model(tiny_hubbard_config(;L,W,boundary_x=bx,boundary_y=by))
        @test length(model.bonds)==expected
        @test all(model.sublattice[i]!=model.sublattice[j] for (i,j) in model.bonds)
    end
    @test_throws ArgumentError hubbard_model(tiny_hubbard_config(L=3,boundary_x=:PBC))
    # Preserve both staggered hopping gauges against the independent Fock matrix.
    A=[0. 0. 1. 0.;1. 0. 0. 0.;0. 1. 0. 0.;0. 0. 0. 1.]
    B=copy(A); B[1,3]=-1
    F=StoMPS.SU2SU2Fermion
    for (pair,first_map,second_map) in ((F.hopping_AB,A,B),(F.hopping_BA,B,A))
        left,right=pair
        parity_right=-right*F.Z
        @tensor hopping[aout bout;ain bin] := left[aout;ain auxiliary]*parity_right[auxiliary bout;bin]
        basis=kron(second_map,first_map)
        @test norm(hopping-hopping')<1e-12
        @test norm(basis*reshape(convert(Array,hopping),16,16)*basis'+dense_two_site_hubbard(interaction=0.))<1e-12
    end
    K=dense_two_site_hubbard(interaction=8.0)
    N=Diagonal([count_ones(i) for i in 0:15])
    H=K+4N
    for symmetry in (:su2su2,:u1su2)
        c=tiny_hubbard_config(;symmetry,ensemble=symmetry==:su2su2 ? :half_filling : :canonical)
        model=hubbard_model(c)
        seed=i->Random.seed!(1200+i)
        thermal=symmetry==:su2su2 ? su2su2_thermal_trajectory(model.H,c;before_step=seed) :
            u1su2_thermal_trajectory(model.H,model.N,c;before_step=seed)
        rho=symmetry==:su2su2 ? thermal.state : thermal.thermal_mpo
        exact_eig=eigvals(Hermitian(K))
        weights=exp.(-c.beta0.*exact_eig)
        if symmetry==:su2su2
            @test thermal.logZ0 ≈ log(sum(weights)) atol=2e-8
            @test thermal.rows[end].energy_per_site ≈ sum(weights.*exact_eig)/sum(weights)/2 atol=2e-8
        else
            @test -thermal.alpha[end] ≈ 4c.beta0 atol=2e-8
            @test thermal.log_normalization[end] ≈ log(sum(weights)) atol=2e-8
            @test thermal.mean_N[end] ≈ 2 atol=2e-9
        end
        X=symmetry==:su2su2 ? vec([collect(x) for x in Iterators.product(1:2,1:2)]) :
            [collect(x) for x in Iterators.product(1:3,1:3) if sum(x.-1)==2]
        born=[hubbard_born(rho,c;rng=MersenneTwister(1),configuration=x) for x in X]
        ps=[hubbard_sample(b,model,c,i) for (i,b) in enumerate(born)]
        @test sum(exp(p.log_p0) for p in ps) ≈ 1 atol=1e-10
        for p in ps
            function step!(state,op,delta,D,index,generation,lineage)
                Random.seed!(1400+lineage)
                env=Environment(state',op,state;disk=false)
                TDVPSweep1!(env,delta;K=32,tol=1e-12,CBEAlg=cbe_algorithm(D,c),
                    trunc=tensor_truncation(D,1e-15),GCstep=false,GCsweep=false,verbose=0)
            end
            advance_hubbard_sample!(p,model,c,2;step!)
        end
        # Enumeration includes p0 explicitly; ordinary Monte Carlo never does.
        ell=[p.log_p0+p.a for p in ps]
        exact=symmetry==:su2su2 ? K : H[findall(==(2),diag(N)),findall(==(2),diag(N))]
        eigenvalues=eigvals(Hermitian(exact)); beta=last(c.beta_list)
        ew=exp.(-beta.*eigenvalues)
        @test logweighted_mean([p.energy for p in ps],ell) ≈ sum(ew.*eigenvalues)/sum(ew) atol=2e-7
        @test logsumexp(ell) ≈ log(sum(ew))-log(sum(exp.(-c.beta0.*eigenvalues))) atol=2e-7
        scores=[configuration_log_probability(p,x,c) for p in ps,x in X]
        @test vec(sum(exp.(scores);dims=2)) ≈ ones(length(ps)) atol=2e-10
        if symmetry==:u1su2
            @test all(p.N==c.Ne for p in ps)
            @test hubbard_expectation(ps[1].state,model.N;disk=false) ≈ c.Ne atol=1e-10
        end
    end
end

@testset "Multiplicity, fugacity and capped inclusion identities" begin
    a,N,m=[1000.,998.,1001.],[0,2,5],[3,1,4]
    fitted=fit_lambda(a,N,m,2.7,-.3;tolerance=1e-12)
    @test fitted.mean_N ≈ 2.7 atol=1e-12
    indices=vcat([fill(i,m[i]) for i in eachindex(m)]...)
    expanded=fit_lambda(a[indices],N[indices],ones(Int,sum(m)),2.7,-.3;tolerance=1e-12)
    @test expanded.lambda ≈ fitted.lambda atol=1e-12
    @test log_weights(sir_weight_base(7.,N,fitted.lambda,-.3),N,fitted.lambda,-.3) ≈ fill(7.,3)
    @test_throws NumberSupportError fit_lambda([0.,0.],[2,3],[1,1],2.,0.)
    @test multiplicity_metrics(a,m,sum(m)).eta ≈ sampling_efficiency(a[indices])
    @test loo_alpha_selection([0.],[0.],[-Inf]).selected_alpha==0
    @test loo_alpha_selection([0.,0.],[-1.,-1.],[-Inf,-Inf]).selected_alpha==0

    masses=[.6,.25,.1,.05]; n=2
    design=sir_inclusion_probabilities(log.(masses),n)
    pi=design.inclusion_probability
    for order in ([1,2,3,4],[4,3,2,1])
        light=[i for i in order if 0<pi[i]<1]
        cuts=sort!(unique(vcat(0.,1.,mod.(cumsum(pi[light]),1.))))
        total=zeros(4); inclusion=zeros(4)
        for j in 1:length(cuts)-1
            probability=cuts[j+1]-cuts[j]
            probability>1e-14 || continue
            selected=sir_systematic_selection(pi,order,(cuts[j]+cuts[j+1])/2)
            @test length(selected)==n && allunique(selected)
            corrected=masses[selected]./pi[selected]
            @test sum(corrected) ≈ 1 atol=1e-14
            inclusion[selected].+=probability
            total[selected].+=probability.*corrected
        end
        @test inclusion ≈ pi atol=1e-14
        @test total ≈ masses atol=1e-14
    end
    plan=sir_pool_draws_without_replacement(log.(masses),n,7;
        order_rng=MersenneTwister(19),offset_rng=MersenneTwister(20))
    ids=plan.active_candidate_id
    @test exp.(plan.candidate_output_log_weight[ids]) ≈ (n/7).*masses[ids]./pi[ids] atol=1e-14
    @test sum(exp.(sir_inclusion_probabilities([0.,-700.,-701.],2).log_inclusion_probability[2:3])) ≈ 1 atol=1e-13

    # An old HT mass is not multiplied by the count of fresh hits. The fresh
    # mass is 2*(p0/Q)*exp(h)=3; the old mass remains 4, giving merged mass 7.
    s=[(;id=1,a=log(4.),m=1,configuration=[1,2],N=2),
       (;id=2,a=0.,m=1,configuration=[2,1],N=2)]
    f=[(;id=3,a=log(.06/.16*4),m=2,configuration=[1,2],N=2)]
    c=(;Ns=2,ensemble=:half_filling,resampling_method=:sir,sir_draw_mode=:without_replacement_configuration)
    selected=select_resampling_pool(s,f,c,(;M=2),0.,0.;rng=MersenneTwister(19))
    @test exp.(selected.audit.mass) ≈ [7.,1.]
    @test exp.(selected.weight) ≈ [3.5,.5]
    @test selected.ids==[1,2] && selected.multiplicity==[1,1]
    @test sir_candidate_count(1.,2)==0
    empty_fresh=fresh_configuration_plan(s,zeros(Int,0,2),merge(c,(;sir_draw_mode=:with_replacement)))
    @test empty_fresh.number_representatives==0 && isempty(empty_fresh.cached_old_id)
    unchanged=select_resampling_pool(s,s[1:0],c,(;M=0),0.,0.;rng=MersenneTwister(19))
    @test exp.(unchanged.weight) ≈ [4.,1.]
end
