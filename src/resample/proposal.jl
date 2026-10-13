"""Draw one member of the empirical mixture for each whole configuration.

Marginalizing this member gives exactly the same mixture law as sequential
posterior-RDM averaging. MPI only has to communicate complete configurations.
"""
function resampling_proposal(s,c,lambda,lambda0,scores; rng::AbstractRNG)
    w = sample_logweights(s,c,lambda,lambda0)
    m = [p.m for p in s]
    Nv = sum(m)
    metric = multiplicity_metrics(w,m,Nv)
    eta_p = clamp(metric.eta_p*(Nv/c.Ns),1/c.Ns,length(s)/c.Ns)
    eta = c.eta_weight == :estimator ? metric.eta :
        sampling_efficiency([p.a-p.log_p0+p.log_Q for p in s])
    mass = w .+ log.(m)
    if c.alpha_selection == :fixed
        alpha = c.resampling_alpha
        selection = (;selected_alpha=alpha)
    else
        loo = [isempty(others) ? -Inf :
            logsumexp(mass[others] .+ scores[others,j])-logsumexp(mass[others])
            for j in eachindex(s) for others in ([i for i in eachindex(s) if i!=j],)]
        amplitude = if c.alpha_selection == :loo_multiplicity
            0.5 .* [log(p.m)+w[i]+p.h+
                (c.ensemble == :grand_canonical ? (lambda-lambda0)*p.N : 0)
                for (i,p) in enumerate(s)]
        else
            [w[i]+0.5*(p.log_Q-p.log_p0) for (i,p) in enumerate(s)]
        end
        selection = loo_alpha_selection(amplitude,[p.log_p0 for p in s],loo;
            alpha_max=c.alpha_max,tolerance=c.alpha_tolerance)
        alpha = selection.selected_alpha
    end
    keep = c.resampling_method == :hybrid ? clamp(floor(Int,eta*c.Ns),1,c.Ns-1) : 0
    M = c.resampling_method == :sir ?
        sir_candidate_count(c.sir_draw_mode == :without_replacement_configuration ?
            metric.eta : eta_p,c.Ns) : c.Ns-keep
    component = [rand(rng) < alpha for _ in 1:M]
    members = [component[i] ? categorical_draw(mass;rng) : 0 for i in 1:M]
    ancestors = [categorical_draw(mass;rng) for _ in 1:keep]
    return (;alpha,selection,M,members,ancestors,eta_pre=eta,
        eta_p_pre=eta_p,common_keep=logsumexp(mass)-log(Nv))
end

"""Prepare fresh unique configurations; reuse old states only in configuration WOR."""
function fresh_configuration_plan(s, X,c)
    fresh = c.resampling_method == :sir && size(X,1)>0 ? deduplicate_configurations(X) :
        (; configurations=copy(X), number_representatives=size(X,1),
            first_draw_id=collect(1:size(X,1)), draw_to_representative=collect(1:size(X,1)),
            multiplicity=ones(Int,size(X,1)))
    old = c.sir_draw_mode == :without_replacement_configuration ?
        [something(findfirst(p->p.configuration==vec(fresh.configurations[i,:]),s),0)
            for i in 1:fresh.number_representatives] :
        zeros(Int,fresh.number_representatives)
    if c.resampling_method == :sir && c.sir_draw_mode == :without_replacement_configuration
        Np = length(union(Set(Tuple(p.configuration) for p in s),Set(Tuple(x) for x in eachrow(X))))
        Np < c.Ns && error("configuration WOR pool has Np=$Np < Ns=$(c.Ns)")
        Np == c.Ns && @warn "configuration WOR retains the entire pool" Np
    end
    return (;fresh...,cached_old_id=old)
end

"""Select a population from scalar pool masses."""
function select_resampling_pool(s,f,c,proposal,lambda,lambda0; rng::AbstractRNG)
    oldw = sample_logweights(s,c,lambda,lambda0)
    freshw = sample_logweights(f,c,lambda,lambda0)
    if c.resampling_method == :regenerate
        ids = [length(s)+i for i in eachindex(f)]
        return (;pool=vcat(s,f),ids,multiplicity=[p.m for p in f],
            weight=freshw,audit=(;kind=:regenerate))
    elseif c.resampling_method == :hybrid
        # Keep each drawn ancestor as a separate lineage copy.
        pool = vcat(s,f)
        ids = vcat(proposal.ancestors,[length(s)+i for i in eachindex(f)])
        return (;pool,ids,multiplicity=vcat(ones(Int,length(proposal.ancestors)),[p.m for p in f]),
            weight=vcat(fill(proposal.common_keep,length(proposal.ancestors)),freshw),
            audit=(;kind=:hybrid,ancestors=proposal.ancestors))
    end
    pool = vcat(s,f)
    mass = vcat(oldw .+ log.([p.m for p in s]),freshw .+ log.([p.m for p in f]))
    if c.sir_draw_mode == :with_replacement_until_unique
        design = _sir_pool_draws_until_unique(vcat(oldw,freshw),[p.m for p in pool],c.Ns;rng)
        ids = design.active_candidate_id
        return (;pool,ids,multiplicity=design.multiplicity[ids],
            weight=design.candidate_output_log_weight[ids],
            audit=(;kind=c.sir_draw_mode,mass,Nv_pre=sum(p.m for p in s),design))
    end
    representatives = collect(eachindex(pool))
    input_to_pool = collect(eachindex(pool))
    if c.sir_draw_mode == :without_replacement_configuration
        grouping = deduplicate_configurations(reduce(vcat,[permutedims(p.configuration) for p in pool]))
        input_to_pool = grouping.draw_to_representative
        representatives = grouping.first_draw_id
        mass = [logsumexp(mass[findall(==(i),input_to_pool)])
            for i in eachindex(representatives)]
    end
    design = if c.sir_draw_mode == :with_replacement
        sir_pool_draws(mass,c.Ns,c.Ns+proposal.M;rng)
    else
        # Independently derived RNGs for the permutation and systematic offset.
        sir_pool_draws_without_replacement(mass,c.Ns,c.Ns+proposal.M;
            order_rng=MersenneTwister(rand(rng,UInt64)),offset_rng=MersenneTwister(rand(rng,UInt64)))
    end
    active = design.active_candidate_id
    weights = c.sir_draw_mode == :with_replacement ?
        fill(design.common_log_weight,length(active)) : design.candidate_output_log_weight[active]
    return (;pool,ids=representatives[active],multiplicity=design.multiplicity[active],
        weight=weights,audit=(;kind=c.sir_draw_mode,mass,input_to_pool,representatives,design))
end
