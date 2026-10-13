# Callers supply before_step(index) to seed FiniteMPS's randomized CBE
# at every thermal node; the shared driver derives these operation seeds.
"""SU2SU2 purification of H-mu*N with explicit application RNG ownership."""
function su2su2_thermal_trajectory(H, c;
        before_step::Function, beta_grid=c.beta_list[c.beta_list .<= c.beta0],
        on_step=(args...)->nothing)
    state = identityMPO(c.L*c.W, su2su2_born_space())
    logZ = normalize_logweight!(state)
    rows = Any[]
    environment = nothing
    try
    for index in eachindex(beta_grid)
        before_step(index)
        started = time()
        if index == 1
            state, _ = SETTN(H, first(beta_grid);
                ρ₀=state, disk=c.environment_disk, CBEAlg=cbe_algorithm(c.settn_bond_dimension,c),
                trunc=tensor_truncation(c.settn_bond_dimension,c.settn_cutoff),
                maxorder=c.settn_maxorder, maxiter=c.settn_maxiter,
                lsnoise=c.settn_noise, tol=c.settn_tolerance,
                compress=c.settn_compression, GCstep=get(c,:settn_gc_step,c.gc_step),
                GCsweep=c.gc_sweep, verbose=get(c,:settn_verbose,0))
            logZ += normalize_logweight!(state)
            diagnostics = (; maximum_truncation=NaN, maximum_entropy=NaN,
                maximum_bond_dimension=state_max_bond_dimension(state))
            elapsed_seconds=time()-started
            environment=Environment(state',H,state;disk=c.environment_disk)
        else
            diagnostics=purification_step!(state,H,purification_method(c),
                beta_grid[index-1],beta_grid[index],c;environment)
            logZ+=diagnostics.log_norm_increment
            elapsed_seconds=diagnostics.elapsed_seconds
        end
        energy = real(scalar!(environment;normalize=true))
        row = (; beta=beta_grid[index], energy_per_site=energy/(c.L*c.W),
            mean_N=Float64(c.L*c.W), mu=c.mu, logZ,
            maximum_truncation=diagnostics.maximum_truncation,
            maximum_entropy=diagnostics.maximum_entropy,
            maximum_bond_dimension=diagnostics.maximum_bond_dimension,
            elapsed_seconds)
        push!(rows,row)
        on_step(index,(; state,rows),row)
    end
    return (; state,rows,beta_index=length(rows),logZ0=logZ,lambda0=last(beta_grid)*c.mu)
    finally
        c.environment_disk && environment!==nothing && _cleanup_environment!(environment)
    end
end

function initial_alpha_mpo(number_sites::Int, target_N::Real)
    isfinite(target_N) && 0 < target_N < 2number_sites || throw(
        ArgumentError("mean particle number must lie strictly between 0 and 2L"))
    target_density = target_N / number_sites
    n_operator = U1SU2Fermion.n
    alpha = log((2-target_density)/target_density)
    rho = MPO(number_sites)
    auxiliary_space = trivial(domain(n_operator)[1])
    for site in 1:number_sites
        rho[site] = permute(
            isometry(auxiliary_space, auxiliary_space) ⊗
                exp(-alpha * n_operator / 2),
            (1, 2),
            (4, 3),
        )
    end
    canonicalize!(rho, 1)
    return rho, alpha
end

function _flow_number!(rho, number_environment, target_N::Real, config)
    number_sites = length(rho)
    tolerance = config.number_tolerance_per_site * number_sites
    canonicalize!(rho, 1)
    canonicalize!(number_environment, 1)
    particle_number = real(scalar!(number_environment; normalize=true))
    abs(particle_number - target_N) <= tolerance && return 0.0, particle_number
    number_gradient = real(connection!(number_environment; moveback=false)[1])
    delta_number = particle_number - target_N
    accumulated_alpha = 0.0
    for _ in 1:config.number_flow_maximum_iterations
        isfinite(number_gradient) && !iszero(number_gradient) || error(
            "particle-number flow has a zero or invalid gradient",
        )
        delta_alpha = delta_number / number_gradient
        exponential = LocalOperator(
            exp(-(delta_alpha / 2) * U1SU2Fermion.n),
            :expn, 0, false, 1,
        )
        sites = Center(rho)[1] == number_sites ?
            reverse(1:number_sites) : (1:number_sites)
        for site in sites
            canonicalize!(rho, site)
            updated = exponential * rho[site]
            local_norm = norm(updated)
            isfinite(local_norm) && local_norm > 0 || error(
                "particle-number flow produced an invalid local norm",
            )
            rmul!(updated, inv(local_norm))
            rmul!(rho, local_norm)
            rho[site] = updated
        end
        accumulated_alpha += delta_alpha
        canonicalize!(number_environment, 1, number_sites)
        final_number = real(scalar!(number_environment; normalize=true))
        number_gradient = (particle_number - final_number) / delta_alpha
        particle_number = final_number
        delta_number = particle_number - target_N
        abs(delta_number) <= tolerance && Center(rho)[2] == 1 && break
    end
    canonicalize!(number_environment, 1, number_sites)
    abs(particle_number - target_N) <= tolerance || error(
        "particle-number flow failed: <N>=$particle_number, target=$target_N",
    )
    return accumulated_alpha, particle_number
end

function _trajectory_state(rho, energy_per_site, mean_number, alpha_history,
        segment_mu, log_normalization_history, diagnostics)
    return (
        thermal_mpo=rho,
        energy_per_site,
        mean_N=mean_number,
        alpha=alpha_history,
        segment_mu,
        log_normalization=log_normalization_history,
        diagnostics,
    )
end

function u1su2_thermal_trajectory(hamiltonian, number_mpo, config;
        before_step::Function, beta_grid=config.beta_list[config.beta_list .<= config.beta0],
        on_step=(args...) -> nothing,
        on_event=(args...) -> nothing)
    energy_per_site = Float64[]
    mean_number = Float64[]
    alpha_history = Float64[]
    segment_mu = Union{Missing,Float64}[]
    log_normalization_history = Float64[]
    diagnostics = Any[]

    on_event(:settn_start, (beta=first(beta_grid),))
    before_step(1)
    settn_started = time()
    rho, alpha_initial = initial_alpha_mpo((config.L*config.W), config.Ne)
    rho, free_energy = SETTN(
        hamiltonian,
        first(beta_grid);
        CBEAlg=cbe_algorithm(config.settn_bond_dimension, config),
        lsnoise=config.settn_noise,
        maxiter=config.settn_maxiter,
        ρ₀=rho,
        disk=config.environment_disk,
        maxorder=config.settn_maxorder,
        verbose=get(config,:settn_verbose,0),
        compress=config.settn_compression,
        GCstep=get(config,:settn_gc_step,config.gc_step),
        GCsweep=config.gc_sweep,
        tol=config.settn_tolerance,
        trunc=tensor_truncation(config.settn_bond_dimension, config.settn_cutoff),
    )
    initial_log_normalization =
        -last(free_energy) * first(beta_grid)
    normalize!(rho)
    environment_H = Environment(rho', hamiltonian, rho;
        disk=config.environment_disk)
    environment_N = Environment(rho', number_mpo, rho;
        disk=config.environment_disk)
    try
        alpha = alpha_initial
        log_normalization = initial_log_normalization
        mean_N = real(scalar!(environment_N; normalize=true))
        if abs(mean_N - config.Ne) >
                config.number_tolerance_per_site * (config.L*config.W)
            correction, mean_N = _flow_number!(rho, environment_N,
                config.Ne, config)
            alpha += correction
            log_normalization += normalize_logweight!(rho)
            Center(environment_H)[:] = [1, length(environment_H)]
            canonicalize!(environment_H, 1)
            canonicalize!(environment_N, 1)
        end
        initialization_seconds=time()-settn_started
        push!(energy_per_site,
            real(scalar!(environment_H; normalize=true)) / (config.L*config.W))
        push!(mean_number, mean_N)
        push!(alpha_history, alpha)
        push!(segment_mu, missing)
        push!(log_normalization_history, log_normalization)
        step_diagnostics = (
            beta=first(beta_grid),
            delta_beta=0.0,
            phase="SETTN",
            maximum_truncation=NaN,
            maximum_entropy=NaN,
            maximum_bond_dimension=state_max_bond_dimension(rho),
            elapsed_seconds=initialization_seconds,
        )
        push!(diagnostics, step_diagnostics)
        on_step(1,
            _trajectory_state(rho, energy_per_site, mean_number,
                alpha_history, segment_mu, log_normalization_history,
                diagnostics),
            step_diagnostics)
        for beta_index in 2:length(beta_grid)
            delta_beta = beta_grid[beta_index] -
                beta_grid[beta_index - 1]
            before_step(beta_index)
            summary=purification_step!(rho,environment_H,environment_N,purification_method(config),
                beta_grid[beta_index-1],beta_grid[beta_index],
                (;energy=energy_per_site[end]*(config.L*config.W),number=mean_number[end],alpha=alpha_history[end]),config)
            log_normalization+=summary.log_norm_increment
            current_N,next_alpha=summary.number,summary.alpha
            push!(energy_per_site,
                real(scalar!(environment_H; normalize=true)) / (config.L*config.W))
            push!(mean_number, current_N)
            push!(segment_mu, -(next_alpha - alpha_history[end]) / delta_beta)
            push!(alpha_history, next_alpha)
            push!(log_normalization_history, log_normalization)
            step_diagnostics = (
                beta=beta_grid[beta_index],
                delta_beta,
                phase="TDVP",
                maximum_truncation=summary.maximum_truncation,
                maximum_entropy=summary.maximum_entropy,
                maximum_bond_dimension=summary.maximum_bond_dimension,
                elapsed_seconds=summary.elapsed_seconds,
            )
            push!(diagnostics, step_diagnostics)
            on_step(beta_index,
                _trajectory_state(rho, energy_per_site, mean_number,
                    alpha_history, segment_mu, log_normalization_history,
                    diagnostics),
                step_diagnostics)
        end
        all(isfinite, energy_per_site) || error("thermal energy is not finite")
        all(isfinite, alpha_history) || error("thermal fugacity is not finite")
        all(isfinite, log_normalization_history) || error(
            "thermal log normalization is not finite",
        )
        all(value -> ismissing(value) || isfinite(value), segment_mu) || error(
            "thermal segment chemical potential is not finite",
        )
        all(value -> abs(value - config.Ne) <=
            config.number_tolerance_per_site * (config.L*config.W), mean_number) ||
            error("thermal trajectory left the configured mean-N tolerance")
        length(diagnostics) == length(beta_grid) || error(
            "thermal diagnostics do not cover the complete beta grid",
        )
        return merge(_trajectory_state(rho, energy_per_site, mean_number,
            alpha_history, segment_mu, log_normalization_history, diagnostics),
            (beta=copy(beta_grid),))
    finally
        config.environment_disk && _cleanup_environment!(environment_H)
        config.environment_disk && _cleanup_environment!(environment_N)
    end
end

function purification_step!(rho::MPO,environment_H,environment_N,::TanTRG,
        beta_from,beta_to,previous,config)
    delta_beta=beta_to-beta_from
    connection=connection!(environment_H,environment_N;moveback=true)
    delta_alpha=-(connection[1,2]*delta_beta+config.Ne-previous.number)/connection[2,2]
    shift=(previous.energy,previous.number)
    started=time_ns()
    info,_=TDVPSweep1!(environment_H,-delta_beta/2,environment_N,-delta_alpha/2;
        K=config.krylov_dimension,tol=config.tdvp_tolerance,
        CBEAlg=cbe_algorithm(config.purification_bond_dimension,config),
        GCstep=config.gc_step,GCsweep=config.gc_sweep,verbose=get(config,:evolution_verbose,0),
        trunc=tensor_truncation(config.purification_bond_dimension,config.purification_cutoff),
        E_shift=shift)
    log_norm_increment=normalize_logweight!(rho)-delta_beta*shift[1]-delta_alpha*shift[2]
    current_N=real(scalar!(environment_N;normalize=true))
    next_alpha=previous.alpha+delta_alpha
    if abs(current_N-config.Ne)>config.number_tolerance_per_site*length(rho)
        correction,current_N=_flow_number!(rho,environment_N,config.Ne,config)
        next_alpha+=correction
        log_norm_increment+=normalize_logweight!(rho)
        Center(environment_H)[:]=[1,length(environment_H)]
        canonicalize!(environment_H,1)
        canonicalize!(environment_N,1)
    end
    elapsed_seconds=(time_ns()-started)/1e9
    return (;tdvp_second_half_diagnostics(info,rho)...,log_norm_increment,
        number=current_N,alpha=next_alpha,elapsed_seconds)
end
