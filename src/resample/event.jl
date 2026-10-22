function score_configurations(particles,s,X,c,e,out)
    score=_phase(e) do
        rows=zeros(length(s),size(X,1))
        for p in particles
            row=findfirst(q->q.id==p.id,s)
            for j in axes(X,1)
                rows[row,j]=configuration_log_probability(p,vec(X[j,:]),c)
                periodic_save!(out;progress=(;phase=:proposal_scoring,sample_index=p.sample_index,configuration=j,ensemble_complete=false))
            end
        end
        rows
    end
    _sum(score,e)
end

function refit(s,c,lambda,lambda0)
    c.ensemble!=:grand_canonical && return lambda
    fit_lambda([p.a for p in s],[p.N for p in s],[p.m for p in s],c.Ne,lambda0;
        initial=lambda,tolerance=c.number_fit_tolerance_per_site*c.L*c.W,
        maximum_iterations=c.number_fit_maximum_iterations).lambda
end

function should_resample(c,stats,beta,events)
    c.resample && length(events)<c.maximum_resampling_events || return false
    c.resampling_policy==:threshold_first && !isempty(events) && return false
    c.resampling_policy==:threshold_until_target && any(e->e.eta_post>=c.eta_target,events) && return false
    c.resampling_policy==:fixed_beta && return beta in c.forced_betas
    c.resampling_policy==:forced_and_threshold && beta in c.forced_betas && return true
    stats.eta<c.eta_threshold
end

function _log_replay!(out,parts,c,lambda,lambda0,logZ0,event)
    out.main===nothing || _section(out.main,"Replay fresh samples from beta0",REPLAY_COLUMNS)
    _note(out,:replay;event,frozen_lambda=lambda,lambda0)
    records=reduce(vcat,parts)
    for beta in sort!(unique(r.beta for r in records))
        rows=filter(r->r.beta==beta,records)
        s=[r.scalar for r in rows];w=sample_logweights(s,c,lambda,lambda0);m=[p.m for p in s]
        errors=conditional_weight_errors(w,m)
        eta=multiplicity_metrics(w,m,sum(m)).eta
        lnZ=logZ0+logsumexp(w.+log.(m))-log(sum(m))
        seconds=maximum(sum(r.seconds for r in rankrows if r.beta==beta;init=0.0) for rankrows in parts)
        row=(;event,beta,lnZ,errors...,eta,elapsed_seconds=seconds)
        push!(out.data["observables"]["replay"],row)
        out.main===nothing || _row(out.main,REPLAY_COLUMNS,(beta,lnZ,errors.lnZ_error,eta,errors.eta_error,seconds))
    end
end

function resample!(particles,s,parent,model,c,e,index,event,lambda,lambda0,logZ0,step!,out)
    started=_now_seconds()
    oldX=reduce(vcat,[permutedims(p.configuration) for p in s])
    scores=c.alpha_selection==:fixed ? zeros(0,0) : score_configurations(particles,s,oldX,c,e,out)
    plan=_broadcast(()->resampling_proposal(s,c,lambda,lambda0,scores;
        rng=rng(c,0x50524f504f53414c,event,0)),e)
    _root(e) do
        _note(out,:Resampling;event,alpha=plan.alpha,loo_seconds=_now_seconds()-started,new_samples=plan.M)
        out.data["diagnostics"]["pending_event"]=(;event,index,plan,pre=s)
    end
    baseline=Dict{Int,Any}()
    X=_phase(e) do
        configurations=zeros(Int,plan.M,c.L*c.W)
        for draw in 1:plan.M
            member=plan.members[draw]
            if member==0 && owns(draw,e)
                born=hubbard_born(parent,c;rng=rng(c,0x4652455348424f52,event,draw))
                baseline[draw]=born;configurations[draw,:]=born.configuration
            elseif member>0
                local_id=findfirst(p->p.id==s[member].id,particles)
                if local_id!==nothing
                    configurations[draw,:]=draw_sample_configuration(particles[local_id],c;
                        rng=rng(c,0x454d504952494341,event,draw))
                end
            end
            periodic_save!(out;progress=(;phase=:proposal_draw,event,draw,ensemble_complete=false))
        end
        configurations
    end
    X=_sum(X,e)
    freshplan=_broadcast(()->fresh_configuration_plan(s,X,c),e)
    uniqueX=freshplan.configurations
    score=score_configurations(particles,s,uniqueX,c,e,out)
    mass=sample_logweights(s,c,lambda,lambda0).+log.([p.m for p in s])
    fresh=HubbardSample[];replay=NamedTuple[]
    base_index=out.sample_counter
    _phase(e) do
        for u in axes(uniqueX,1)
            old=freshplan.cached_old_id[u];draw=freshplan.first_draw_id[u]
            cached=old==0 ? nothing : findfirst(p->p.id==s[old].id,particles)
            (old==0 ? owns(draw,e) : cached!==nothing) || continue
            log_Pr=logsumexp(mass.+score[:,u])-logsumexp(mass)
            if old==0
                born=haskey(baseline,draw) ? baseline[draw] :
                    hubbard_born(parent,c;rng=rng(c,0x524547454e455241,event,draw),configuration=vec(uniqueX[u,:]))
                Q=linear_mixture_log_probability(born.log_probability,log_Pr,plan.alpha)
                p=hubbard_sample(born,model,c,length(s)+u;generation=event,
                    lineage=event*(c.Ns+1)+draw,m=freshplan.multiplicity[u],log_Q=Q)
                p.sample_index=base_index+u
                for k in 1:index
                    d=k==1 ? initial_diagnostics(p.state) : advance_hubbard_sample!(p,model,c,k;step!)
                    beta=sample_betas(c)[k]
                    _sample_log(out,p,beta,d,only(sample_logweights([sample_scalars(p)],c,lambda,lambda0));phase=:replay)
                    push!(replay,(;beta,scalar=sample_scalars(p),seconds=d.elapsed_seconds))
                    _particle_progress!(out,p,beta,:replay)
                end
            else
                ancestor=particles[cached]
                Q=linear_mixture_log_probability(ancestor.log_p0,log_Pr,plan.alpha)
                p=HubbardSample(;id=length(s)+u,sample_index=base_index+u,
                    lineage=event*(c.Ns+1)+draw,generation=event,state=ancestor.state,
                    configuration=copy(ancestor.configuration),N=ancestor.N,energy=ancestor.energy,
                    h=ancestor.h,m=freshplan.multiplicity[u],a=ancestor.log_p0-Q+ancestor.h,
                    log_p0=ancestor.log_p0,log_Q=Q,compression_norm=ancestor.compression_norm,
                    diagnostics=ancestor.diagnostics)
            end
            push!(fresh,p)
        end
    end
    fs=global_scalars(fresh,e)
    replay_parts=_gather(replay,e)
    _root(e) do
        _log_draws(out,[(;index=draw,log_probability=fs[freshplan.draw_to_representative[draw]].log_Q)
            for draw in 1:plan.M])
        _log_replay!(out,replay_parts,c,lambda,lambda0,logZ0,event)
    end
    selection=_broadcast(()->select_resampling_pool(s,fs,c,plan,lambda,lambda0;
        rng=rng(c,0x53495253454c4543,event,0)),e)
    localpool=vcat(particles,fresh);next=HubbardSample[]
    _phase(e) do
        for (newid,pool_index) in enumerate(selection.ids)
            record=selection.pool[pool_index]
            local_index=findfirst(p->p.id==record.id,localpool)
            local_index===nothing && continue
            p=deepcopy(localpool[local_index]);p.id=newid;p.m=selection.multiplicity[newid]
            p.sample_index=base_index+length(fs)+newid
            p.a=selection.weight[newid]-(c.ensemble==:grand_canonical ? (lambda-lambda0)*p.N : 0)
            push!(next,p)
        end
    end
    out.sample_counter=base_index+length(fs)+length(selection.ids)
    post=global_scalars(next,e)
    lambda_post=_broadcast(e) do
        validate_population(post,c)
        refit(post,c,lambda,lambda0)
    end
    eta_post=_broadcast(()->sample_statistics(post,c,lambda_post,lambda0).eta,e)
    accepted=c.acceptance_policy==:always || eta_post>plan.eta_pre
    audit=(;index,beta=sample_betas(c)[index],event,plan,freshplan,configurations=X,
        scores,proposal_scores=score,selection=selection.audit,pre=s,fresh=fs,post,
        lambda_pre=lambda,lambda_post,eta_post,accepted,
        evolved_fresh_representatives=count(iszero,freshplan.cached_old_id),
        reused_old_configurations=count(!iszero,freshplan.cached_old_id))
    _root(e) do
        delete!(out.data["diagnostics"],"pending_event")
        _note(out,:resampling_result;event,accepted,eta_post)
        out.main===nothing || _section(out.main,"stochastic MPS",STO_COLUMNS)
    end
    accepted ? (next,post,lambda_post,audit) : (particles,s,lambda,audit)
end
