using Random

mutable struct PoolDrawTape <: AbstractRNG
    values::Vector{Float64}
    position::Int
end
PoolDrawTape(values) = PoolDrawTape(Float64.(values),0)
function Random.rand(rng::PoolDrawTape)
    rng.position += 1
    return rng.values[rng.position]
end

@testset "Adaptive pool draws retain multiplicities and mean-weight scale" begin
    draw = StoMPS._sir_pool_draws_until_unique
    # Equal masses, unequal per-copy weights and input multiplicities.
    w,m = log.([1.,2.,4.]),[4,2,1]
    rng=PoolDrawTape([.1,.1,.4,.9])
    d=draw(w,m,3;rng)
    @test d.selected_candidate_id == [1,1,2,3]
    @test d.multiplicity == [2,1,1]
    @test (d.draw_count,d.Nv_post,d.stop_reason) == (4,4,:target_reached)
    @test exp.(d.candidate_output_log_weight) ≈ fill(12/7,3)
    @test rng.position == 4
    capped=draw(w,m,3;rng=PoolDrawTape(fill(.1,9)))
    @test capped.multiplicity == [9,0,0]
    @test (capped.draw_count,capped.stop_reason) == (9,:draw_limit)
    lastdraw=draw(zeros(2),ones(Int,2),2;rng=PoolDrawTape([0.,0.,0.,0.,0.,.9]))
    @test lastdraw.stop_reason == :target_reached && lastdraw.draw_count == 6

    # Keeping a small pool is deterministic and need not satisfy the draw cap.
    rng=PoolDrawTape(Float64[])
    kept=draw(log.([2.,5.]),[8,4],3;rng)
    @test kept.active_candidate_id == [1,2] && kept.multiplicity == [8,4]
    @test kept.candidate_output_log_weight == log.([2.,5.])
    @test kept.Nv_post == kept.Nv_pool == 12
    @test kept.stop_reason == :pool_under_target && kept.draw_count == rng.position == 0

    s=[(;id=1,a=log(2.),m=3,configuration=[1,2],N=2,energy=-2.,h=.2,log_p0=-1.,log_Q=-1.),
       (;id=2,a=log(1.),m=2,configuration=[2,1],N=2,energy=1.,h=.1,log_p0=-1.,log_Q=-1.)]
    f=[merge(s[1],(;id=3,m=1,a=log(3.)))]
    c=(;Ns=3,ensemble=:canonical,resampling_method=:sir,
        sir_draw_mode=:with_replacement_until_unique,eta_weight=:estimator,
        alpha_selection=:loo_multiplicity,alpha_max=.999,alpha_tolerance=1e-8)
    selected=select_resampling_pool(s,f,c,(;M=1),0.,0.;rng=PoolDrawTape([.1,.6,.9]))
    @test selected.ids == [1,2,3] # IDs 1 and 3 have the same configuration.
    @test selected.audit.Nv_pre == 5 && selected.audit.design.Nv_pool == 6
    @test exp.(selected.weight) ≈ fill(11/6,3)
    small=select_resampling_pool(s,s[1:0],c,(;M=0),0.,0.;rng=PoolDrawTape(Float64[]))
    @test small.multiplicity == [3,2] && small.weight == [p.a for p in s]

    stats=sample_statistics(s,c,0.,0.)
    expanded=[1,1,1,2,2]
    reference=multiplicity_weighted_statistics([p.energy for p in s][expanded],
        [p.a for p in s][expanded],ones(Int,5),5)
    @test stats.Nv == 5 && stats.K == 2 && stats.Np == 2
    @test stats.energy ≈ reference.mean && stats.stderr ≈ reference.errorbar
    @test stats.eta ≈ reference.effective_sample_size/5
    @test stats.eta_p ≈ 8^2/(3*(6^2+2^2))
    proposal=resampling_proposal(s,c,0.,0.,fill(-log(2.),2,2);rng=MersenneTwister(33))
    @test proposal.M == 3-floor(Int,8^2/(6^2+2^2))
    @test proposal.eta_p_pre ≈ stats.eta_p
    @test logsumexp(log.([p.m for p in s]).+[p.a for p in s])-log(stats.Nv) ≈ log(8/5)

    # Legacy fixed-count sampling is still the same sequence of iid draws.
    legacy=sir_pool_draws(log.([.8,.2]),2,6;rng=PoolDrawTape([.1,.1]))
    @test legacy.selected_candidate_id == [1,1] && legacy.multiplicity == [2,0]
    @test legacy.common_log_weight ≈ -log(6)
end

@testset "Exact stopped two-ID distribution" begin
    p=.8
    observed=0.; normalization=0.
    # Enumerate fixed six-draw paths independently of when the implementation stops.
    for sequence in Iterators.product(fill(1:2,6)...)
        probability=prod(i==1 ? p : 1-p for i in sequence)
        tape=PoolDrawTape([i==1 ? p/2 : (1+p)/2 for i in sequence])
        d=StoMPS._sir_pool_draws_until_unique(log.([p,1-p]),ones(Int,2),2;rng=tape)
        observed += probability*d.multiplicity[1]/d.Nv_post
        normalization += probability
    end
    # The first change at n has mass p^(n-1)(1-p) or (1-p)^(n-1)p.
    exact=p^6+sum(p^(n-1)*(1-p)*(n-1)/n+(1-p)^(n-1)*p/n for n in 2:6)
    @test normalization ≈ 1
    @test observed ≈ exact ≈ .716992
    @test abs(observed-p) > .08
end
