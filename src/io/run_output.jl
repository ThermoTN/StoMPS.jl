mutable struct RunOutput
    config
    execution::AbstractExecution
    run_id::String
    data::Dict{String,Any}
    main::Union{Nothing,IOStream}
    sample_ids::Set{Int}
    sample_phase::Dict{Int,Symbol}
    sample_counter::Int
    result_path::Union{Nothing,String}
    last_save::Float64
end
_now_seconds()=time_ns()/1e9
_column_width(name)=max(length(name),name=="probability" ? 24 : 20)
function _table_header(io,columns)
    header=join((lpad(name,_column_width(name)) for name in columns)," | ")
    println(io,header);println(io,repeat('-',length(header)));flush(io)
end
function _section(io,title,columns)
    println(io,"======================== ",title," ================================")
    _table_header(io,columns)
end
function _parameters(io,values,prefix="")
    for key in sort!(collect(keys(values));by=string)
        value=values[key];name=isempty(prefix) ? string(key) : prefix*"."*string(key)
        if value isa NamedTuple || value isa AbstractDict
            _parameters(io,value,name)
        else
            println(io,'"',name,"\": ",repr(value))
        end
    end
end
function _row(io,columns,values)
    cells=(v isa AbstractString ? v : v isa Integer ? string(v) :
        @sprintf("%.10e",v===nothing || ismissing(v) ? NaN : Float64(v)) for v in values)
    println(io,join((lpad(value,_column_width(name)) for (name,value) in zip(columns,cells))," | "))
    flush(io)
end
const PUR_COLUMNS=("beta","energy_total","energy_per_site","lnZ","max_MPO_truncation","time(s)")
const STO_COLUMNS=("beta","energy_total","energy_per_site","energy_per_site_error","lnZ","lnZ_error","eta","eta_error","max_MPS_truncation","time(s)")
const SAMPLE_COLUMNS=("beta","energy_total","energy_per_site","log_weight","max_MPS_truncation","count","time(s)")
const SAMPLE_REPLAY_COLUMNS=("beta","log_weight","max_MPS_truncation","count","time(s)")
const REPLAY_COLUMNS=("beta","lnZ","lnZ_error","eta","eta_error","time(s)")
const DRAW_COLUMNS=("index","probability")

function _open_output(c,e,stage,model;checkpoint=nothing,run_id=nothing)
    id=if run_id===nothing
        _broadcast(e) do
            if checkpoint!==nothing && c.output_dir!==nothing
                path=joinpath(c.output_dir,"results","result.jld2")
                if isfile(path)
                    old=load_result(path)["metadata"]
                    old["stage"]=="purification" && old["stage_complete"] && old["run_id"]==checkpoint.run_id && return checkpoint.run_id
                end
            end
            Dates.format(now(UTC),dateformat"yyyymmddTHHMMSS")*"-"*string(time_ns())
        end
    else
        run_id
    end
    data=Dict{String,Any}("parameters"=>c,"observables"=>Dict{String,Any}(
        "curve"=>NamedTuple[],"purification"=>NamedTuple[],"replay"=>NamedTuple[]),
        "diagnostics"=>Dict{String,Any}("events"=>Any[]),
        "metadata"=>Dict{String,Any}("run_id"=>id,"stage"=>String(stage),
            "complete"=>false,"stage_complete"=>false,"restart_supported"=>false,
            "mpi_size"=>_nprocs(e),"julia_threads"=>Threads.nthreads(),
            "blas_threads"=>LinearAlgebra.BLAS.get_num_threads(),"progress"=>nothing,
            "uncertainty"=>"paired delta-method estimates for lnZ and eta"))
    checkpoint===nothing || append!(data["observables"]["purification"],checkpoint.rows)
    checkpoint===nothing || (data["metadata"]["checkpoint_run_id"]=checkpoint.run_id)
    checkpoint===nothing || (data["metadata"]["checkpoint_provenance"]=checkpoint.provenance)
    main=nothing;path=nothing
    _root(e) do
        data["metadata"]["provenance"]=_source_record(c,model)
        if c.output_dir!==nothing
            stem=stage in (:purification,:stochastic) ? "result" : String(stage)*"_result"
            path=joinpath(c.output_dir,"results",stem*".jld2")
            logpath=joinpath(c.output_dir,"log","run.log")
            append_log=false
            if isfile(path)
                old=load_result(path)
                append_log=stage==:stochastic && old["metadata"]["stage"]=="purification" &&
                    old["metadata"]["stage_complete"] && old["metadata"]["run_id"]==id
                append_log || c.overwrite_output || error("output exists; choose a new output_dir or overwrite_output=true")
            elseif isfile(logpath)
                append_log=stage in (:dynamics,:ground_state,:ed,:purification_reference)
                append_log || c.overwrite_output || error("run.log exists without a matching completed preparation")
            end
            mkpath(dirname(logpath));mkpath(dirname(path))
            main=open(logpath,append_log ? "a" : "w")
            println(main,"============================== Parameters ==================================")
            _parameters(main,c)
            _parameters(main,(;stage,run_id=id,mpi_size=_nprocs(e),julia_threads=Threads.nthreads(),
                blas_threads=LinearAlgebra.BLAS.get_num_threads(),
                julia_version=string(VERSION),stomps_version=string(Base.pkgversion(StoMPS)),
                finite_mps_version=string(Base.pkgversion(FiniteMPS))))
            println(main,"# lnZ/eta errors use paired delta-method estimates")
            println(main,repeat('=',76));flush(main)
        end
    end
    out=RunOutput(c,e,id,data,main,Set{Int}(),Dict{Int,Symbol}(),c.Ns,path,_now_seconds())
    _root(()->periodic_save!(out;force=true),e)
    return out
end

function periodic_save!(out;force=false,progress=nothing)
    _rank(out.execution)==0 || return
    progress===nothing || (out.data["metadata"]["progress"]=progress)
    out.result_path===nothing && return
    current=_now_seconds()
    force || current-out.last_save>=out.config.save_interval_seconds || return
    d=out.data
    save_result(out.result_path;parameters=d["parameters"],observables=d["observables"],
        diagnostics=d["diagnostics"],metadata=d["metadata"])
    out.last_save=current
end

function _finish_output!(out;complete=false,stage_complete=false,error=nothing)
    if _rank(out.execution)==0
        out.data["metadata"]["complete"]=complete
        out.data["metadata"]["stage_complete"]=stage_complete
        error===nothing || (out.data["metadata"]["error"]=error)
        try
            periodic_save!(out;force=true)
        finally
            out.main===nothing || close(out.main)
        end
    end
end

function _note(out,stage;fields...)
    out.main===nothing && return
    print(out.main,"# ",stage)
    for (key,value) in pairs(fields)
        print(out.main,' ',key,'=',repr(value))
    end
    println(out.main);flush(out.main)
end

function _probability_string(logp)
    logp==-Inf && return "0.0000000000e+00"
    exponent=floor(Int,logp/log(10.0))
    @sprintf("%.10fe%+d",exp(logp-exponent*log(10.0)),exponent)
end
function _log_draws(out,draws;header=true)
    out.main===nothing && return
    header && _table_header(out.main,DRAW_COLUMNS)
    for r in sort(draws;by=r->r.index)
        _row(out.main,DRAW_COLUMNS,(r.index,_probability_string(r.log_probability)))
    end
end
function _log_purification(out,r)
    out.main===nothing && return
    _row(out.main,PUR_COLUMNS,(r.beta,r.energy_per_site*out.config.L*out.config.W,
        r.energy_per_site,r.logZ,r.maximum_truncation,r.elapsed_seconds))
    _note(out,:purification;mean_N=r.mean_N,mu=r.mu,D=r.maximum_bond_dimension,entropy=r.maximum_entropy,
        wall_elapsed_seconds=r.wall_elapsed_seconds,cumulative_evolution_seconds=r.cumulative_evolution_seconds)
end
function _log_ensemble(out,r;phase=:ordinary)
    out.main===nothing && return
    phase==:pre_resampling && println(out.main,"# before resampling")
    phase==:post_resampling && println(out.main,"# after resampling; errors condition on inherited common scale")
    _row(out.main,STO_COLUMNS,(r.beta,r.energy_per_site*out.config.L*out.config.W,
        r.energy_per_site,r.stderr_per_site,r.logZ,r.lnZ_error,r.eta,r.eta_error,
        r.maximum_truncation,r.elapsed_seconds))
    _note(out,:population;Ns=r.Ns,Nv=r.Nv,K=r.K,Np=r.Np,eta_p=r.eta_p,
        mean_N=r.mean_N,mu=r.mu,D=r.maximum_bond_dimension,entropy=r.maximum_entropy)
end

function _sample_log(out,p,beta,diagnostics,log_weight;phase=:ordinary)
    out.config.output_dir===nothing && return
    directory=joinpath(out.config.output_dir,"log","samples");mkpath(directory)
    path=joinpath(directory,@sprintf("rank_%03d_sample_%03d.log",_rank(out.execution),p.sample_index))
    first=!(p.sample_index in out.sample_ids)
    columns=phase==:replay ? SAMPLE_REPLAY_COLUMNS : SAMPLE_COLUMNS
    open(path,first ? "w" : "a") do io
        if first
            println(io,"# run_id=",out.run_id," seed=",out.config.master_seed," sample_index=",p.sample_index,
                " lineage=",p.lineage," generation=",p.generation," configuration=",repr(p.configuration))
        end
        if first || get(out.sample_phase,p.sample_index,:none)!=phase
            println(io,"# phase=",phase)
            _section(io,"finite-MPS",columns)
        end
        values=phase==:replay ? (beta,log_weight,diagnostics.maximum_truncation,p.m,diagnostics.elapsed_seconds) :
            (beta,p.energy,p.energy/(out.config.L*out.config.W),log_weight,diagnostics.maximum_truncation,p.m,diagnostics.elapsed_seconds)
        _row(io,columns,values)
        println(io,"# physical_log_norm=",p.h," D=",diagnostics.maximum_bond_dimension," entropy=",diagnostics.maximum_entropy)
    end
    push!(out.sample_ids,p.sample_index);out.sample_phase[p.sample_index]=phase
end

function _particle_progress!(out,p,beta,phase)
    if _rank(out.execution)==0
        out.data["diagnostics"]["in_progress_sample"]=(;phase,beta,scalar=sample_scalars(p))
    end
    periodic_save!(out;progress=(;phase,beta,sample_index=p.sample_index,ensemble_complete=false))
end
