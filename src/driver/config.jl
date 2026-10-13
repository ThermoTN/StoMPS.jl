const RUN_KEYS=(:output_dir,:overwrite_output,:reuse_beta0,:save_interval_seconds,
    :stop_after_beta,:task,:expected_mpi_size,:expected_julia_threads,:blas_threads)
identity_config(c;exclude=RUN_KEYS)=NamedTuple{Tuple(k for k in keys(c) if k ∉ exclude)}(
    Tuple(c[k] for k in keys(c) if k ∉ exclude))
dynamics_source_identity(c)=identity_config(c;exclude=(RUN_KEYS...,
    :dynamics_betas,:dynamics_dt,:dynamics_tmax,:frequency_list,:window_cutoff,:broadening,
    :ground_state_dynamics,:ground_state_twice_spin,:ground_state_bond_dimensions,
    :ground_state_sweeps,:ground_state_energy_tolerance))

"""Resolve optional application controls once, before logging or execution."""
function effective_config(c)
    haskey(c,:resume) && error("resume was removed in 0.6; use reuse_beta0 and restart from beta0")
    merge((;lattice_kind=:square,purification=TanTRG(),evolution=TDVP(),
        output_dir=nothing,overwrite_output=false,reuse_beta0=true,save_interval_seconds=3600.0,
        stop_after_beta=nothing,task=:all,expected_mpi_size=nothing,expected_julia_threads=nothing,
        blas_threads=1,resampling_alpha=.5,resampling_policy=:threshold_each,
        forced_betas=Float64[],acceptance_policy=:always,eta_target=.9,
        measure_local_observables=false,dynamics_betas=Float64[],ground_state_dynamics=false,
        number_fit_tolerance_per_site=1e-8,number_fit_maximum_iterations=128),c)
end

function validate_config(c,e::AbstractExecution=SerialExecution())
    all(isfinite,c.beta_list) && all(>(0),c.beta_list) && all(diff(c.beta_list).>0) ||
        error("beta_list must be strictly increasing and positive")
    c.beta0 in c.beta_list && c.beta0>first(c.beta_list) || error("beta0 must be a purification endpoint")
    c.Ns>1 && 0<c.Ne<2c.L*c.W || error("invalid sample count or target particle number")
    c.ensemble==:canonical && !(c.Ne isa Integer) && error("canonical Ne must be an integer")
    all(>(0),(c.settn_bond_dimension,c.purification_bond_dimension,c.sample_bond_dimension)) || error("bond dimensions must be positive")
    0<=c.eta_threshold<=1 && 0<=c.resampling_alpha<=1 && 0<c.alpha_max<1 && c.alpha_tolerance>0 || error("invalid resampling controls")
    c.sir_draw_mode in (:with_replacement,:with_replacement_until_unique,:without_replacement_configuration,:without_replacement_pool_id) || error("unknown SIR draw mode")
    c.resampling_method in (:sir,:hybrid,:regenerate) || error("unknown resampling method")
    c.alpha_selection in (:fixed,:loo_first,:loo_generalized,:loo_multiplicity) || error("unknown alpha selector")
    c.resampling_policy in (:threshold_each,:threshold_first,:fixed_beta,:threshold_until_target,:forced_and_threshold) || error("unknown resampling policy")
    c.acceptance_policy in (:always,:eta_improvement) || error("unknown acceptance policy")
    c.eta_weight in (:estimator,:physical_norm) || error("unknown eta weights")
    if c.symmetry!=:su2su2
        c.resampling_method==:sir && c.sir_draw_mode in (:with_replacement,:with_replacement_until_unique) &&
            c.alpha_selection==:loo_multiplicity && c.eta_weight==:estimator &&
            c.acceptance_policy==:always && c.resampling_policy==:threshold_each ||
            error("U1SU2 supports with-replacement multiplicity LOO-SIR")
    end
    c.resampling_method==:sir && (c.alpha_selection!=:loo_multiplicity || c.eta_weight!=:estimator) && error("SIR requires estimator eta and multiplicity LOO")
    c.resampling_method!=:sir && c.sir_draw_mode!=:with_replacement && error("this draw mode requires SIR")
    c.resampling_method==:hybrid && c.alpha_selection!=:fixed && error("hybrid requires fixed alpha")
    c.resampling_method==:regenerate && c.alpha_selection==:loo_multiplicity && error("regeneration requires fixed or nonmultiplicity LOO")
    c.alpha_selection==:loo_first && c.maximum_resampling_events!=1 && error("first-event LOO requires one event")
    c.acceptance_policy==:eta_improvement && !(c.resampling_method==:regenerate && c.alpha_selection==:loo_generalized) && error("eta acceptance requires generalized-LOO regeneration")
    c.maximum_resampling_events>=0 || error("negative event budget")
    all(b->b in sample_betas(c),vcat(c.dynamics_betas,c.forced_betas)) || error("requested betas must belong to the sample grid")
    if !isempty(c.dynamics_betas) || c.task in (:dynamics,:ed) || c.ground_state_dynamics
        c.symmetry==:u1su2 && c.ensemble==:canonical && c.W==1 && c.boundary_x==c.boundary_y==:OBC || error("JP dynamics requires a canonical U1SU2 open chain")
    end
    c.task==:ground_state && c.ensemble!=:canonical && error("DMRG requires canonical U1SU2")
    c.save_interval_seconds>0 || error("save interval must be positive")
    c.expected_mpi_size===nothing || c.expected_mpi_size==_nprocs(e) || error("MPI size mismatch")
    c.expected_julia_threads===nothing || c.expected_julia_threads==Threads.nthreads() || error("Julia thread count mismatch")
    c.purification isa PurificationMethod && c.evolution isa EvolutionMethod || error("invalid evolution method")
    return c
end

function _prepare_run(c,comm)
    e=_execution(comm)
    c=_phase(()->validate_config(effective_config(c),e),e)
    signatures=_gather(repr(c),e)
    _root(()->all(==(first(signatures)),signatures) || error("MPI ranks supplied different configurations"),e)
    LinearAlgebra.BLAS.set_num_threads(c.blas_threads)
    return c,e
end
