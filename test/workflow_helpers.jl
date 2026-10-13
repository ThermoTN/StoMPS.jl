function workflow_config(output;application="U1SU2_Hubbard_square_lattice_canonical_ensemble",kwargs...)
    defaults=Module(gensym(:Defaults))
    Base.include(defaults,joinpath(@__DIR__,"..","examples",application,"config.jl"))
    base=Base.invokelatest(getfield(defaults,:config))
    merge(base,(;L=2,W=1,Ne=2,Ns=32,boundary_x=:OBC,boundary_y=:OBC,
        beta_list=[1/128,1/16,1/8,1/4],beta0=1/16,
        settn_bond_dimension=32,purification_bond_dimension=32,sample_bond_dimension=32,
        eta_threshold=1.0,alpha_max=.9,output_dir=output,environment_disk=true,
        gc_step=false,gc_sweep=false,task=:all,save_interval_seconds=1e-6,
        dynamics_betas=Float64[],ground_state_bond_dimensions=[32],ground_state_sweeps=12),
        (;kwargs...))
end

function contains_tensor(x)
    x isa StoMPS.FiniteMPS.AbstractMPS && return true
    x isa StoMPS.FiniteMPS.TensorKit.AbstractTensorMap && return true
    x isa AbstractDict && return any(contains_tensor,values(x))
    x isa NamedTuple || x isa Tuple || x isa AbstractArray ? any(contains_tensor,x) : false
end

function check_workflow(result,c)
    rows=result["observables"]["curve"]
    @test result["metadata"]["complete"]
    @test length(rows)==3
    @test length(result["diagnostics"]["events"])==2
    @test !contains_tensor(result)
    @test all(r->isfinite(r.logZ) && isfinite(r.energy_per_site) && sum(p.m for p in r.scalars)==r.Nv,rows)
    @test all(r->r.K<=c.Ns && r.Np==length(Set(Tuple(p.configuration) for p in r.scalars)),rows)
    c.ensemble==:canonical && @test all(p.N==c.Ne for r in rows for p in r.scalars)
    for row in rows
        expanded=reduce(vcat,[fill(w,p.m) for (w,p) in zip(row.log_weight,row.scalars)])
        @test row.logZ ≈ estimate_logpartition(result["metadata"]["logZ0"],expanded)
    end
    log=read(joinpath(c.output_dir,"log","run.log"),String)
    @test occursin("# before resampling",log) && occursin("# after resampling",log)
    for (directory,_,files) in walkdir(c.output_dir), file in files
        endswith(file,".jld2") || continue
        path=joinpath(directory,file)
        file=="checkpoint.jld2" && continue
        occursin(joinpath("dynamics","sources"),path) && continue
        file=="ground_state.jld2" && continue
        @test !contains_tensor(load_result(path))
    end
    @test !isfile(joinpath(c.output_dir,"resume.jld2"))
end

function compare_curves(reference,other)
    @test length(reference)==length(other)
    for (a,b) in zip(reference,other)
        @test (a.beta,a.Ns,a.Nv,a.K,a.Np)==(b.beta,b.Ns,b.Nv,b.K,b.Np)
        @test [a.energy_per_site,a.logZ,a.eta,a.eta_p,a.mean_N] ≈
            [b.energy_per_site,b.logZ,b.eta,b.eta_p,b.mean_N] atol=2e-9 rtol=2e-9
        @test [(p.configuration,p.m,p.lineage,p.generation) for p in a.scalars]==
            [(p.configuration,p.m,p.lineage,p.generation) for p in b.scalars]
        @test [p.a for p in a.scalars] ≈ [p.a for p in b.scalars] atol=2e-9 rtol=2e-9
        @test [p.h for p in a.scalars] ≈ [p.h for p in b.scalars] atol=2e-9 rtol=2e-9
    end
end
