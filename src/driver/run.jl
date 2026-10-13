initial_diagnostics(state)=(;maximum_truncation=NaN,maximum_entropy=NaN,
    maximum_bond_dimension=state_max_bond_dimension(state),elapsed_seconds=0.0)

function validate_population(s,c)
    !isempty(s) && all(p->p.m isa Integer && p.m>0,s) || error("invalid population multiplicities")
    Nv=sum(p.m for p in s)
    if c.sir_draw_mode==:with_replacement_until_unique
        length(s)<=c.Ns && Nv>=c.Ns || error("invalid adaptive population")
    else
        Nv==c.Ns || error("multiplicities do not sum to Ns")
    end
    all(p->isfinite(p.energy) && isfinite(p.a) && isfinite(p.h),s) || error("nonfinite sample")
    c.ensemble==:canonical && !all(p->p.N==c.Ne,s) && error("canonical particle number mismatch")
    Nv
end

function _sample_step!(state,H,delta,D,index,generation,lineage,c)
    seed=deterministic_seed(c.master_seed,0x544456505f434245,generation,lineage*1_000_000+index)
    evolution_step!(state,H,evolution_method(c),delta,c;rng=MersenneTwister(seed),seed,D,
        cutoff=c.sample_cutoff)
end

"""Prepare beta0 once on rank zero. MPI callers load the published MPO locally."""
function run_purification(config;comm=nothing,reference=false,before_step=nothing)
    c,e=_prepare_run(config,comm)
    path=c.output_dir===nothing ? nothing : joinpath(c.output_dir,"results","checkpoint.jld2")
    existing=_broadcast(()->path!==nothing && isfile(path),e)
    if existing
        c.reuse_beta0 || c.overwrite_output || error("beta0 checkpoint exists; reuse it or explicitly overwrite output")
        if c.reuse_beta0
            cached=_phase(()->load_checkpoint(path,c),e)
            !reference && return cached
        end
    end
    c.output_dir===nothing && _nprocs(e)>1 && error("MPI preparation requires a shared output_dir")
    model=_phase(()->hubbard_model(c;lattice_kind=c.lattice_kind),e)
    out=_open_output(c,e,reference ? :purification_reference : :purification,model)
    cp=nothing
    try
        _root(e) do
            trajectory_started=_now_seconds()
            cumulative_evolution_seconds=0.0
            out.main===nothing || _section(out.main,"finite-MPO",PUR_COLUMNS)
            grid=reference ? c.beta_list : c.beta_list[c.beta_list.<=c.beta0]
            rows=out.data["observables"]["purification"]
            seed_step=before_step===nothing ?
                i->Random.seed!(deterministic_seed(c.thermal_seed,0x544845524d414c43,0,i)) : before_step
            function on_step(i,s,d)
                row=c.symmetry==:su2su2 ? d :
                    (;beta=grid[i],energy_per_site=s.energy_per_site[end],mean_N=s.mean_N[end],
                        mu=-s.alpha[end]/grid[i],logZ=s.log_normalization[end],
                        maximum_truncation=d.maximum_truncation,maximum_entropy=d.maximum_entropy,
                        maximum_bond_dimension=d.maximum_bond_dimension,elapsed_seconds=d.elapsed_seconds)
                cumulative_evolution_seconds+=row.elapsed_seconds
                row=merge(row,(;wall_elapsed_seconds=_now_seconds()-trajectory_started,
                    cumulative_evolution_seconds))
                push!(rows,row);_log_purification(out,row)
                periodic_save!(out;progress=(;phase=:purification,beta=grid[i],ensemble_complete=true))
                if grid[i]==c.beta0
                    state=c.symmetry==:su2su2 ? s.state : s.thermal_mpo
                    cp=ThermalCheckpoint(state,row.logZ,row.mu*c.beta0,thermal_identity(c),
                        copy(rows),out.run_id,out.data["metadata"]["provenance"])
                    path===nothing || (existing && c.reuse_beta0) || save_checkpoint(path,cp)
                end
            end
            thermal_trajectory(purification_method(c),model,c;before_step=seed_step,beta_grid=grid,on_step)
        end
        _root(()->_finish_output!(out;complete=reference,stage_complete=true),e)
    catch exception
        _finish_output!(out;error=sprint(showerror,exception))
        rethrow()
    end
    reference && return _rank(e)==0 ? out.data : nothing
    path===nothing ? cp : _phase(()->load_checkpoint(path,c),e)
end

"""A purification method may override initialization and the full endpoint schedule."""
function thermal_trajectory(::TanTRG,model,c;kwargs...)
    c.symmetry==:su2su2 ? su2su2_thermal_trajectory(model.H,c;kwargs...) :
        u1su2_thermal_trajectory(model.H,model.N,c;kwargs...)
end

function _ensemble_diagnostics(local_diagnostics,e)
    records=_gather(local_diagnostics,e)
    _broadcast(e) do
        rows=reduce(vcat,records)
        # NaN is an undefined diagnostic, not a numeric initial maximum.
        function measured_max(key)
            values=[getproperty(r,key) for r in rows if isfinite(getproperty(r,key))]
            isempty(values) ? NaN : maximum(values)
        end
        (;maximum_truncation=measured_max(:maximum_truncation),
            maximum_entropy=measured_max(:maximum_entropy),
            maximum_bond_dimension=isempty(rows) ? 0 : maximum(r.maximum_bond_dimension for r in rows),
            elapsed_seconds=maximum(sum(r.elapsed_seconds for r in local_rows;init=0.0) for local_rows in records))
    end
end

function _ensemble_row(s,c,lambda,lambda0,logZ0,beta,d;local_observables=nothing,correlations=nothing)
    validate_population(s,c)
    stats=sample_statistics(s,c,lambda,lambda0)
    logZ=logZ0+logsumexp(stats.w.+log.(stats.m))-log(stats.Nv)
    isfinite(logZ) || error("nonfinite logZ")
    (;beta,energy_per_site=stats.energy/(c.L*c.W),stderr_per_site=stats.stderr/(c.L*c.W),
        mean_N=stats.mean_N,mu=c.symmetry==:su2su2 ? c.mu : c.ensemble==:grand_canonical ? lambda/beta : nothing,
        eta=stats.eta,eta_p=stats.eta_p,Ns=c.Ns,Nv=stats.Nv,K=stats.K,Np=stats.Np,logZ,
        lnZ_error=stats.lnZ_error,eta_error=stats.eta_error,lambda,
        log_weight=copy(stats.w),eta_weight=c.eta_weight,
        physical_log_norm=[p.h for p in s],scalars=s,local_observables,correlations,d...)
end

function local_measurements(particles,s,c,e,stats)
    c.measure_local_observables || return nothing
    rows=_phase(e) do
        values=zeros(length(s),2c.L*c.W)
        for p in particles
            v=hubbard_local_observables(p.state;disk=c.environment_disk)
            values[findfirst(q->q.id==p.id,s),:]=vcat(v.number,v.double)
        end
        values
    end
    rows=_sum(rows,e)
    _broadcast(e) do
        values=[c.ensemble==:grand_canonical ?
            fitted_statistics(rows[:,i],[p.N for p in s],stats.w,stats.m) :
            multiplicity_weighted_statistics(rows[:,i],stats.w,stats.m,stats.Nv) for i in axes(rows,2)]
        (;number=[v.mean for v in values[1:c.L*c.W]],double=[v.mean for v in values[c.L*c.W+1:end]],
            number_stderr=[v.errorbar for v in values[1:c.L*c.W]],double_stderr=[v.errorbar for v in values[c.L*c.W+1:end]])
    end
end

"""Run Born sampling and cooling from beta0, returning the numerical history on rank zero."""
function run_stochastic(config;checkpoint=nothing,checkpoint_path=nothing,comm=nothing)
    c,e=_prepare_run(config,comm)
    model=_phase(()->hubbard_model(c;lattice_kind=c.lattice_kind),e)
    cp=_phase(e) do
        if checkpoint===nothing
            path=checkpoint_path===nothing ? (c.output_dir===nothing ? nothing :
                joinpath(c.output_dir,"results","checkpoint.jld2")) : checkpoint_path
            path===nothing && error("supply checkpoint or checkpoint_path")
            load_checkpoint(path,c)
        else
            checkpoint.preparation==thermal_identity(c) || error("incompatible thermal checkpoint")
            isapprox(norm(checkpoint.state),1;atol=c.probability_tolerance,rtol=0) || error("thermal checkpoint is not normalized")
            checkpoint
        end
    end
    out=_open_output(c,e,:stochastic,model;checkpoint=cp)
    particles=HubbardSample[];lambda=cp.lambda0;lambda0=cp.lambda0;logZ0=cp.logZ0
    history=out.data["observables"]["curve"];events=out.data["diagnostics"]["events"]
    step! = (state,H,d,D,i,g,l)->_sample_step!(state,H,d,D,i,g,l,c)
    try
        measurement=_phase(e) do
            if c.lattice_kind==:triangular
                bulk=bulk_geometry(c,model);schedule=correlation_schedule(c)
                (;bulk,schedule,operators=isempty(schedule) ? nothing : bulk_operators(c,model,bulk))
            else
                nothing
            end
        end
        sector=_phase(e) do
            sector_weight=0.0
            for id in 1:c.Ns
                owns(id,e) || continue
                born=hubbard_born(cp.state,c;rng=rng(c,0x494e4954424f524e,0,id))
                push!(particles,hubbard_sample(born,model,c,id))
                id==1 && c.ensemble==:canonical && (sector_weight=born.target_sector_weight)
                _particle_progress!(out,particles[end],c.beta0,:born)
            end
            [sector_weight]
        end
        sector=_sum(sector,e)
        c.ensemble==:canonical && (logZ0+=log(sector[1])-lambda0*c.Ne)
        s=global_scalars(particles,e)
        _root(e) do
            out.data["metadata"]["logZ0"]=logZ0
            out.data["metadata"]["lambda0"]=lambda0
            measurement===nothing || (out.data["metadata"]["measurement"]=(;measurement.bulk,measurement.schedule,definition=CORRELATION_DEFINITION))
            out.main===nothing || _section(out.main,"Born sampling",DRAW_COLUMNS)
            _log_draws(out,[(;index=p.sample_index,log_probability=p.log_p0) for p in s];header=false)
            out.main===nothing || _section(out.main,"stochastic MPS",STO_COLUMNS)
        end
        betas=sample_betas(c)
        for index in eachindex(betas)
            beta=betas[index]
            local_diagnostics=_phase(e) do
                diagnostics=Dict{Int,Any}()
                for p in particles
                    d=index==1 ? initial_diagnostics(p.state) : advance_hubbard_sample!(p,model,c,index;step!)
                    diagnostics[p.sample_index]=d
                    _particle_progress!(out,p,beta,:stochastic)
                end
                diagnostics
            end
            s=global_scalars(particles,e)
            lambda=_broadcast(()->refit(s,c,lambda,lambda0),e)
            stats=_broadcast(()->sample_statistics(s,c,lambda,lambda0),e)
            d=_ensemble_diagnostics(collect(values(local_diagnostics)),e)
            trigger=index>1 && should_resample(c,stats,beta,events)
            _phase(e) do
                for p in particles
                    _sample_log(out,p,beta,local_diagnostics[p.sample_index],
                        only(sample_logweights([sample_scalars(p)],c,lambda,lambda0)))
                end
            end
            if trigger
                pre=_broadcast(()->_ensemble_row(s,c,lambda,lambda0,logZ0,beta,d),e)
                _root(()->_log_ensemble(out,pre;phase=:pre_resampling),e)
                particles,s,lambda,audit=resample!(particles,s,cp.state,model,c,e,index,length(events)+1,lambda,lambda0,logZ0,step!,out)
                push!(events,audit)
                stats=_broadcast(()->sample_statistics(s,c,lambda,lambda0),e)
                post_d=_ensemble_diagnostics([p.diagnostics===nothing ? initial_diagnostics(p.state) : p.diagnostics for p in particles],e)
                d=merge(post_d,(;elapsed_seconds=d.elapsed_seconds))
                _phase(e) do
                    for p in particles
                        _sample_log(out,p,beta,initial_diagnostics(p.state),
                            only(sample_logweights([sample_scalars(p)],c,lambda,lambda0));phase=:post_resampling)
                    end
                end
            end
            local_values=local_measurements(particles,s,c,e,stats)
            correlations=measurement!==nothing && any(r->r.index==index,measurement.schedule) ?
                correlation_measurements(particles,s,c,e,stats,measurement.operators,measurement.bulk,index,out) : nothing
            row=_broadcast(()->_ensemble_row(s,c,lambda,lambda0,logZ0,beta,d;
                local_observables=local_values,correlations),e)
            push!(history,row)
            _root(e) do
                delete!(out.data["diagnostics"],"in_progress_sample")
                _log_ensemble(out,row;phase=trigger ? :post_resampling : :ordinary)
                correlations===nothing || _note(out,:correlation;beta,css_tot=correlations.css_tot,stderr=correlations.css_tot_stderr)
                periodic_save!(out;progress=(;phase=:stochastic,beta,ensemble_complete=true))
            end
            beta in c.dynamics_betas && save_dynamics_source!(particles,s,c,e,index,lambda,lambda0,out.run_id)
            c.stop_after_beta!==nothing && beta>=c.stop_after_beta && break
        end
        complete=length(history)==length(sample_betas(c))
        _root(()->_finish_output!(out;complete,stage_complete=complete),e)
        return _rank(e)==0 ? out.data : nothing
    catch exception
        _finish_output!(out;error=sprint(showerror,exception))
        rethrow()
    end
end

"""Dispatch application tasks to the shared thermal, stochastic and reference workflows."""
function run_application(config;comm=nothing)
    c=effective_config(config)
    c.task==:ground_state && return run_ground_state(c;comm)
    c.task==:dynamics && return run_dynamics(c;comm)
    c.task==:ed && return run_ed(c;comm)
    c.task==:purification_reference && return run_purification(c;comm,reference=true)
    c.task in (:all,:thermal,:cooling) || error("unknown task: $(c.task)")
    cp=run_purification(c;comm)
    c.task==:thermal && return cp
    result=run_stochastic(c;checkpoint=cp,comm)
    if c.task==:all && !isempty(c.dynamics_betas) &&
            (c.stop_after_beta===nothing || c.stop_after_beta>=last(sample_betas(c)))
        run_dynamics(c;comm)
    end
    result
end
