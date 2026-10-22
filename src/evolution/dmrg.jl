using FiniteMPS.TensorKit

function initial_ground_state(c)
    n=c.L*c.W
    occupancies=zeros(Int,n)
    initial_rng=rng(c,0x47524f554e445354,0,0)
    for _ in 1:c.Ne
        candidates=findall(==(minimum(occupancies)),occupancies)
        occupancies[rand(initial_rng,candidates)]+=1
    end
    charges=reverse(cumsum(reverse(1 .- occupancies)))
    spin=c.ground_state_twice_spin//2
    left=Rep[U₁ × SU₂]((n-c.Ne,spin)=>1)
    bonds=[i==1 ? left : Rep[U₁ × SU₂]((q,s)=>1 for s in 0:1//2:max(1//2,spin))
        for (i,q) in enumerate(charges)]
    Random.seed!(deterministic_seed(c.master_seed,0x47524f554e445354,0,0))
    return randMPS(Float64,repeat([U1SU2Fermion.pspace],n),bonds)
end

"""Two-site reference in one fixed (Ne,S) sector; convergence refers to sweeps at the final D."""
function dmrg_ground_state(model,c;on_step=(args...)->nothing)
    mod(c.Ne,2)==mod(c.ground_state_twice_spin,2) || error("incompatible Ne and total spin")
    state=initial_ground_state(c)
    history=NamedTuple[]
    schedule=repeat(c.ground_state_bond_dimensions;inner=c.ground_state_sweeps)
    converged=false
    for (sweep,D) in enumerate(schedule)
    # Seed the DMRG initialization and sweep streams separately.
        c.lattice_kind==:triangular && Random.seed!(deterministic_seed(c.master_seed,0x4753444d52475357,0,sweep))
        started=time_ns()
        info=with_environment(state,model.H;disk=c.environment_disk) do env
            first(DMRGSweep2!(env;K=c.krylov_dimension,tol=c.tdvp_tolerance,
                trunc=tensor_truncation(D,c.sample_cutoff),noise=0.0,
                GCstep=c.gc_step,GCsweep=c.gc_sweep,verbose=0))
        end
        canonicalize!(state,1);normalize!(state)
        elapsed_seconds=(time_ns()-started)/1e9
        E=hubbard_expectation(state,model.H;disk=c.environment_disk)
        row=(;sweep,D,energy_per_site=E/(c.L*c.W),
            maximum_truncation=maximum(x.Bond.TrunErr^2 for half in info for x in half),
            maximum_entropy=maximum(x.Bond.SE for half in info for x in half),elapsed_seconds)
        push!(history,row);on_step(row)
        finalrows=filter(r->r.D==last(c.ground_state_bond_dimensions),history)
        converged=length(finalrows)>=4 && all(abs(finalrows[j].energy_per_site-finalrows[j-1].energy_per_site)<c.ground_state_energy_tolerance
            for j in length(finalrows)-2:length(finalrows))
        converged && break
    end
    converged || error("DMRG sweep limit reached; numerical diagnostics retained")
    N=hubbard_expectation(state,model.N;disk=c.environment_disk)
    abs(N-c.Ne)<1e-7 || error("ground-state particle-number mismatch")
    values=hubbard_local_observables(state;disk=c.environment_disk)
    summary=(;energy_per_site=history[end].energy_per_site,particle_number=N,
        twice_spin=c.ground_state_twice_spin,local_observables=values,
        sweep_converged=true,bond_dimension_convergence_verified=false,spin_sector_convergence_verified=false)
    (;state,history,summary)
end
