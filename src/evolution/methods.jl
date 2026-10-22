"""Dispatch boundary for preparing a normalized thermal MPO at requested beta endpoints."""
abstract type PurificationMethod end
"""SETTN initialization followed by one-site tanTRG evolution."""
struct TanTRG <: PurificationMethod end

"""Dispatch boundary for finite-MPS evolution; the step is an exponential time."""
abstract type EvolutionMethod end
"""FiniteMPS one-site TDVP with controlled bond expansion."""
struct TDVP <: EvolutionMethod end

purification_method(c) = get(c, :purification, TanTRG())
evolution_method(c) = get(c, :evolution, TDVP())
cbe_dimension(D, c) = D + floor(Int, D*c.cbe_ratio)
cbe_algorithm(D, c) = NaiveCBE(cbe_dimension(D,c), c.cbe_tolerance; rsvd=c.rsvd)
tensor_truncation(D, cutoff) = truncdim(D) & truncbelow(cutoff)

function with_environment(f, state::DenseMPS, operator; disk=false)
    environment = Environment(state', operator, state; disk)
    try
        return f(environment)
    finally
        disk && _cleanup_environment!(environment)
    end
end

function _cleanup_environment!(environment)
    cleanup!(environment.El)
    cleanup!(environment.Er)
    return nothing
end

"""Normalize a state and return its squared-norm logarithm, counted exactly once."""
function normalize_logweight!(state::DenseMPS)
    n = norm(state)
    isfinite(n) && n > 0 || error("evolution produced a zero or nonfinite norm")
    increment = 2log(n)
    normalize!(state)
    return increment
end

"""
    evolution_step!(state, H, method, delta, config; rng, seed=nothing, D, cutoff)

Evolve and normalize a `DenseMPS`. Return `log_norm_increment` and numerical
diagnostics. The caller accumulates squared-norm and estimator log weights. For thermal samples
pass `delta=-delta_beta/2`. Before calling randomized FiniteMPS operations,
seed Julia's default RNG with `seed`, or draw a seed from `rng` when omitted.
The shared driver derives operation seeds independently of MPI rank.
"""
function evolution_step!(state::DenseMPS, H, ::TDVP, delta, c;
        rng::AbstractRNG, seed=nothing, D=c.sample_bond_dimension, cutoff=c.sample_cutoff)
    Random.seed!(seed===nothing ? rand(rng, UInt64) : seed)
    return _tdvp_step!(state,H,delta,c;D,cutoff)
end

    # Seed this operation from the driver seed or the supplied RNG.
function _tdvp_step!(state::DenseMPS,H,delta,c;D,cutoff)
    with_environment(state,H;disk=c.environment_disk) do env
        _tdvp_environment_step!(state,env,delta,c;D,cutoff)
    end
end

function _tdvp_environment_step!(state::DenseMPS,env,delta,c;D,cutoff)
    started=time_ns()
    info,_=TDVPSweep1!(env,delta;K=c.krylov_dimension,tol=c.tdvp_tolerance,
        CBEAlg=cbe_algorithm(D,c),trunc=tensor_truncation(D,cutoff),
        GCstep=c.gc_step,GCsweep=c.gc_sweep,verbose=get(c,:evolution_verbose,0))
    elapsed_seconds=(time_ns()-started)/1e9
    diagnostics=tdvp_second_half_diagnostics(info,state)
    log_norm_increment=normalize_logweight!(state)
    return (;diagnostics...,log_norm_increment,elapsed_seconds)
end

"""
    purification_step!(state, H, method, beta_from, beta_to, config)

Endpoint-based purification interface. A new method may use its own internal
temperature schedule. Return the squared-norm log increment and diagnostics;
leave the thermal factor normalized. U1SU2 tanTRG has a specialized method
taking the energy and particle-number environments and its thermodynamic state.
"""
function purification_step!(state::MPO,H,::TanTRG,beta_from,beta_to,c;environment=nothing)
    controls=(;D=c.purification_bond_dimension,cutoff=c.purification_cutoff)
    environment===nothing ? _tdvp_step!(state,H,-(beta_to-beta_from)/2,c;controls...) :
        _tdvp_environment_step!(state,environment,-(beta_to-beta_from)/2,c;controls...)
end
