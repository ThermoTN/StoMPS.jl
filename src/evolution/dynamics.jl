"""Total JP greater correlation. Ground states use the connected P source."""
conductivity_particle(p,c,model,snapshot_index;kwargs...)=
    conductivity_particle(evolution_method(c),p,c,model,snapshot_index;kwargs...)

# This real-time pair retains physical amplitudes and the H-E phase convention.
function conductivity_particle(::TDVP,p,c,model,snapshot_index;ground_state=false,on_step=(args...)->nothing)
    dt,tmax=c.dynamics_dt,c.dynamics_tmax
    nsteps=round(Int,tmax/dt)
    dt>0 && tmax>=0 && isapprox(nsteps*dt,tmax;atol=1e-12,rtol=0) || error("tmax must be a nonnegative multiple of dt")
    times=collect(0:nsteps).*dt
    u=complex(deepcopy(p.state)); v=deepcopy(u)
    positions=polarization_positions(model.lattice)
    bonds=current_bonds(model.lattice,positions;t=c.t)
    P=polarization_mpo(model.lattice,positions)
    J=current_mpo(model.lattice,bonds)
    cbe=cbe_algorithm(c.sample_bond_dimension,c)
    trunc=tensor_truncation(c.sample_bond_dimension,c.sample_cutoff)
    function seed!(domain,index)
        Random.seed!(deterministic_seed(c.master_seed,domain,p.generation,
            p.lineage*1_000_000_000+snapshot_index*1_000_000+index))
    end
    seed!(0x44594e414d494350,0)
    mul!(v,P,u;CBEAlg=cbe,trunc,tol=c.tdvp_tolerance,GCsweep=c.gc_sweep,verbose=0)
    polarization_mean=ground_state ? hubbard_expectation(u,P;disk=c.environment_disk) : 0.0
    ground_state && axpy!(-polarization_mean,u,v)
    observable=current_observable_tree(model.lattice,bonds)
    function observe()
        # Preserve the FiniteMPS mixed-bra coefficient convention and serial
        # observable contraction (its threaded LRU cache can miss spuriously).
        calObs!(observable,u,v;serial=true,disk=false)
        correction=conj(coef(u))/coef(u)
        return sum(-im*dx*correction*(t*observable.Refs["fwd_$(i)_$(j)"][(i,j)][]-
            conj(t)*observable.Refs["rev_$(i)_$(j)"][(j,i)][]) for (i,j,dx,t) in bonds)
    end
    correlation=ComplexF64[observe()]
    direct=deepcopy(v); seed!(0x44594e414d49434a,0)
    mul!(direct,J,v;CBEAlg=cbe,trunc,tol=c.tdvp_tolerance,GCsweep=c.gc_sweep,verbose=0)
    t0_direct=inner(u,direct); t0_residual=correlation[1]-t0_direct
    envu=Environment(u',model.H,u;disk=c.environment_disk)
    envv=Environment(v',model.H,v;disk=c.environment_disk)
    E=hubbard_expectation(u,model.H;disk=c.environment_disk)
    diagnostics=Any[]
    try
        for index in 1:nsteps
            started=time()
            seed!(0x44594e414d494355,index)
            du=if ground_state
                nothing
            else
                info,_=TDVPSweep1!(envu,-im*dt;K=c.krylov_dimension,tol=c.tdvp_tolerance,
                    CBEAlg=cbe,trunc,E_shift=E,GCstep=c.gc_step,GCsweep=c.gc_sweep,verbose=0)
                tdvp_second_half_diagnostics(info,u)
            end
            seed!(0x44594e414d494356,index)
            info,_=TDVPSweep1!(envv,-im*dt;K=c.krylov_dimension,tol=c.tdvp_tolerance,
                CBEAlg=cbe,trunc,E_shift=E,GCstep=c.gc_step,GCsweep=c.gc_sweep,verbose=0)
            dv=tdvp_second_half_diagnostics(info,v)
            seconds=time()-started
            push!(correlation,observe())
            push!(diagnostics,(;time=times[index+1],u=du,v=dv,seconds))
            on_step((;time=times[index+1],correlation=correlation[end],diagnostics=diagnostics[end]),
                (;time=times[1:index+1],Cgreater_JP=copy(correlation)))
        end
    finally
        c.environment_disk && (_cleanup_environment!(envu);_cleanup_environment!(envv))
    end
    all(isfinite,correlation) || error("nonfinite JP correlation")
    return (;time=times,Cgreater_JP=correlation,t0_direct,t0_residual,
        t0_residual_abs=abs(t0_residual),diagnostics,positions,bonds,
        response_kind=:greater_JP,response_normalization=:total_extensive,
        P_convention="sum x_j*n_j with centered x",J_convention="+i[H,P]",
        energy_shift=E,ground_state,polarization_mean)
end
