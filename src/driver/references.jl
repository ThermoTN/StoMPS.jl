function run_ground_state(config;comm=nothing)
    c,e=_prepare_run(merge(config,(;task=:ground_state)),comm)
    model=_phase(()->hubbard_model(c;lattice_kind=c.lattice_kind),e)
    out=_open_output(c,e,:ground_state,model)
    try
        _root(e) do
            rows=NamedTuple[];out.data["observables"]["sweeps"]=rows
            columns=("sweep","D","energy_per_site","max_MPS_truncation","entropy","time(s)")
            out.main===nothing || _section(out.main,"DMRG",columns)
            result=dmrg_ground_state(model,c;on_step=r->begin
                push!(rows,r)
                out.main===nothing || _row(out.main,columns,Tuple(r))
                periodic_save!(out;progress=(;phase=:dmrg,sweep=r.sweep))
            end)
            out.data["observables"]["ground_state"]=result.summary
            if c.output_dir!==nothing
                atomic_save(joinpath(c.output_dir,"results","ground_state.jld2");
                    format=RESULT_FORMAT,complete=true,config=identity_config(c),
                    state=result.state,summary=result.summary,run_id=out.run_id)
            end
            if c.ground_state_dynamics
                p=(;state=result.state,id=1,sample_index=1,lineage=1,generation=0)
                out.data["observables"]["dynamics"]=conductivity_particle(p,c,model,0;
                    ground_state=true,on_step=(r,history)->begin
                        _log_dynamics_sample(out,p,r,0)
                        out.data["diagnostics"]["in_progress_dynamics"]=history
                        periodic_save!(out;progress=(;phase=:ground_state_dynamics,time=r.time))
                    end)
                delete!(out.data["diagnostics"],"in_progress_dynamics")
            end
        end
        _root(()->_finish_output!(out;complete=true,stage_complete=true),e)
    catch exception
        _finish_output!(out;error=sprint(showerror,exception));rethrow()
    end
    _rank(e)==0 ? out.data : nothing
end

function _log_dynamics_sample(out,p,r,index)
    out.config.output_dir===nothing && return
    directory=joinpath(out.config.output_dir,"log","samples");mkpath(directory)
    path=joinpath(directory,@sprintf("rank_%03d_sample_%03d.log",_rank(out.execution),p.sample_index))
    first=!(p.sample_index in out.sample_ids)
    columns=("time","real_JP","imag_JP","max_MPS_truncation","D","time(s)")
    # Dynamics appends to the matching stochastic sample history. A fresh
    # standalone run starts its own section without erasing that history.
    open(path,"a") do io
        if first || get(out.sample_phase,p.sample_index,:none)!=Symbol("dynamics_",index)
            println(io,"# run_id=",out.run_id," phase=dynamics snapshot_index=",index,
                " sample_index=",p.sample_index," lineage=",p.lineage)
            _section(io,"JP dynamics",columns)
        end
        _row(io,columns,(r.time,real(r.correlation),imag(r.correlation),
            r.diagnostics.v.maximum_truncation,r.diagnostics.v.maximum_bond_dimension,r.diagnostics.seconds))
    end
    push!(out.sample_ids,p.sample_index);out.sample_phase[p.sample_index]=Symbol("dynamics_",index)
end

"""Compute JP responses from only the explicitly requested thermal MPS snapshots."""
function run_dynamics(config;comm=nothing)
    c,e=_prepare_run(merge(config,(;task=:dynamics)),comm)
    c.output_dir===nothing && error("dynamics requires saved source snapshots")
    model=_phase(()->hubbard_model(c;lattice_kind=c.lattice_kind),e)
    out=_open_output(c,e,:dynamics,model)
    try
        for beta in c.dynamics_betas
            index=findfirst(==(beta),sample_betas(c))
            manifest=joinpath(c.output_dir,"dynamics","sources","beta_$(index)","manifest.jld2")
            source=_broadcast(e) do
                d=load_result(manifest)
                d["complete"] && d["config"]==dynamics_source_identity(c) && d["beta"]==beta || error("incompatible dynamics source")
                d
            end
            directory=joinpath(dirname(manifest),source["directory"])
            local_results=_phase(e) do
                records=NamedTuple[]
                for (position,s) in enumerate(source["scalars"])
                    owns(position,e) || continue
                    saved=load_result(joinpath(directory,"sample_$(s.sample_index).jld2"))
                    saved["record"]==s || error("dynamics sample disagrees with source manifest")
                    p=(;s...,state=saved["state"])
                    result=conductivity_particle(p,c,model,index;on_step=(r,history)->begin
                        _log_dynamics_sample(out,p,r,index)
                        if _rank(e)==0
                            out.data["diagnostics"]["in_progress_dynamics"]=(;beta,sample_index=s.sample_index,history)
                            periodic_save!(out;progress=(;phase=:dynamics,beta,time=r.time,ensemble_complete=false))
                        end
                    end)
                    push!(records,(;id=s.id,sample_index=s.sample_index,result))
                end
                records
            end
            parts=_gather(local_results,e)
            _root(e) do
                records=sort!(reduce(vcat,parts);by=r->r.id)
                s=source["scalars"]
                w=sample_logweights(s,c,source["lambda"],source["lambda0"]);m=[p.m for p in s]
                correlations=reduce(hcat,[r.result.Cgreater_JP for r in records])
                mass=exp.(w.+log.(m).-logsumexp(w.+log.(m)))
                weighted=correlations*mass;times=records[1].result.time;Nv=sum(m)
                transform=isempty(c.frequency_list) ? nothing : OpticalED.finite_time_conductivity(
                    times,weighted,c.L,c.frequency_list;cutoff=c.window_cutoff,eta=c.broadening)
                stderr_real=[multiplicity_weighted_statistics(real.(correlations[i,:]),w,m,Nv).errorbar for i in axes(correlations,1)]
                stderr_imag=[multiplicity_weighted_statistics(imag.(correlations[i,:]),w,m,Nv).errorbar for i in axes(correlations,1)]
                result=(;beta,time=times,Cgreater_JP=weighted,stderr_real,stderr_imag,transform,
                    Ns=c.Ns,Nv,K=length(m),Np=length(Set(Tuple(p.configuration) for p in s)),
                    t0_residual_abs=maximum(r.result.t0_residual_abs for r in records),
                    source_run_id=source["run_id"],samples=records)
                out.data["observables"]["beta_$(index)"]=result
                delete!(out.data["diagnostics"],"in_progress_dynamics")
                columns=("time","real_JP","imag_JP","stderr_real","stderr_imag")
                out.main===nothing || _section(out.main,"JP dynamics beta=$(beta)",columns)
                for i in eachindex(times)
                    out.main===nothing || _row(out.main,columns,(times[i],real(weighted[i]),imag(weighted[i]),stderr_real[i],stderr_imag[i]))
                end
                periodic_save!(out;progress=(;phase=:dynamics,beta,ensemble_complete=true))
            end
        end
        _root(()->_finish_output!(out;complete=true,stage_complete=true),e)
    catch exception
        _finish_output!(out;error=sprint(showerror,exception));rethrow()
    end
    _rank(e)==0 ? out.data : nothing
end

function run_ed(config;comm=nothing)
    c,e=_prepare_run(merge(config,(;task=:ed)),comm)
    c.L<=8 || error("full-block ED is restricted to L<=8")
    model=_phase(()->hubbard_model(c;lattice_kind=c.lattice_kind),e)
    out=_open_output(c,e,:ed,model)
    try
        _root(e) do
            out.data["observables"]["reference"]=OpticalED.block_canonical_optical_ed(c.L,c.Ne;
                hopping=c.t,interaction=c.U,beta=isempty(c.dynamics_betas) ? [c.beta0] : c.dynamics_betas,
                time=collect(0:c.dynamics_dt:c.dynamics_tmax),frequencies=c.frequency_list,
                cutoff=c.window_cutoff,eta=c.broadening,on_block=(phase,data)->begin
                    _note(out,:ED;phase,data)
                    periodic_save!(out;progress=(;phase=:ed,block=data.block_index))
                end)
        end
        _root(()->_finish_output!(out;complete=true,stage_complete=true),e)
    catch exception
        _finish_output!(out;error=sprint(showerror,exception));rethrow()
    end
    _rank(e)==0 ? out.data : nothing
end
