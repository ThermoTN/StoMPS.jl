const RESULT_FORMAT="stomps-0.6-v1"

function atomic_save(path;fields...)
    mkpath(dirname(path))
    temporary=path*".tmp"
    JLD2.jldsave(temporary;fields...)
    mv(temporary,path;force=true)
    return path
end

"""Save numerical parameters, observables, diagnostics and metadata atomically."""
function save_result(path;parameters,observables,diagnostics,metadata)
    atomic_save(path;format=RESULT_FORMAT,parameters,observables,diagnostics,metadata)
end

function load_result(path)
    data=JLD2.load(path)
    get(data,"format",nothing)==RESULT_FORMAT || error("unsupported StoMPS result format")
    return data
end

const THERMAL_KEYS=(:symmetry,:L,:W,:boundary_x,:boundary_y,:lattice_kind,:t,:U,:Ne,
    :settn_bond_dimension,:purification_bond_dimension,:cbe_ratio,:cbe_tolerance,:rsvd,
    :settn_maxorder,:settn_maxiter,:settn_noise,:settn_tolerance,:settn_compression,
    :settn_cutoff,:purification_cutoff,:krylov_dimension,:tdvp_tolerance,
    :number_tolerance_per_site,:number_flow_maximum_iterations,:thermal_seed,:purification)
function thermal_identity(c)
    pairs=Pair{Symbol,Any}[k=>c[k] for k in THERMAL_KEYS]
    c.symmetry==:su2su2 && push!(pairs,:mu=>c.mu)
    (;pairs...,beta0=c.beta0,beta_grid=c.beta_list[c.beta_list.<=c.beta0],
        settn_gc_step=get(c,:settn_gc_step,c.gc_step),
        thermal_seed_schedule=get(c,:thermal_seed_schedule,:derived))
end

struct ThermalCheckpoint{S<:MPO,P}
    state::S
    logZ0::Float64
    lambda0::Float64
    preparation::NamedTuple
    rows::Vector{NamedTuple}
    run_id::String
    provenance::P
end

function save_checkpoint(path,cp::ThermalCheckpoint)
    atomic_save(path;format=RESULT_FORMAT,kind=:beta0,complete=true,state=cp.state,
        logZ0=cp.logZ0,lambda0=cp.lambda0,preparation=cp.preparation,rows=cp.rows,
        run_id=cp.run_id,provenance=cp.provenance)
end

function load_checkpoint(path,c)
    c=effective_config(c)
    d=load_result(path)
    d["kind"]==:beta0 && d["complete"] && d["preparation"]==thermal_identity(c) ||
        error("incompatible beta0 checkpoint")
    d["state"] isa MPO && isapprox(norm(d["state"]),1;atol=c.probability_tolerance,rtol=0) ||
        error("checkpoint must contain a normalized thermal MPO")
    ThermalCheckpoint(d["state"],d["logZ0"],d["lambda0"],d["preparation"],
        NamedTuple[d["rows"]...],d["run_id"],d["provenance"])
end

function _source_record(c,model)
    root=dirname(@__DIR__)
    paths=sort([joinpath(directory,name) for base in (root,joinpath(root,"..","ext"))
        for (directory,_,names) in walkdir(base) for name in names if endswith(name,".jl")])
    digest=bytes2hex(SHA.sha256(join([relpath(p,root)*"\n"*read(p,String) for p in paths],"\n")))
    (;source_digest=digest,julia_version=string(VERSION),
        stomps_version=string(Base.pkgversion(StoMPS)),finite_mps_version=string(Base.pkgversion(FiniteMPS)),
        bonds=model.bonds,positions=model.positions,
        sector_order=c.symmetry==:su2su2 ? string.(collect(sectors(su2su2_born_space()))) : ["empty","single","double"],
        energy_convention=c.symmetry==:su2su2 ? "H-mu*N; no constant shift" : "H only")
end

function save_dynamics_source!(particles,scalars,c,e,index,lambda,lambda0,run_id)
    c.output_dir===nothing && error("dynamics source snapshots require output_dir")
    directory=joinpath(c.output_dir,"dynamics","sources","beta_$(index)",run_id)
    _phase(e) do
        for p in particles
            atomic_save(joinpath(directory,"sample_$(p.sample_index).jld2");
                format=RESULT_FORMAT,state=p.state,record=sample_scalars(p))
        end
    end
    _root(e) do
        all(isfile(joinpath(directory,"sample_$(p.sample_index).jld2")) for p in scalars) || error("incomplete dynamics source")
        atomic_save(joinpath(dirname(directory),"manifest.jld2");format=RESULT_FORMAT,
            complete=true,config=dynamics_source_identity(c),beta=sample_betas(c)[index],index,
            run_id,directory=basename(directory),scalars,lambda,lambda0)
    end
end
