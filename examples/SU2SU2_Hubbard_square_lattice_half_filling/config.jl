# All run parameters live here. Copy to ignored config.local.jl for personal runs.
# From this directory: mpiexecjl --project=. -n 4 julia --project=. -t 2 main.jl
using StoMPS: TanTRG, TDVP

function config()
    L, W = 8, 4
    t, U = 1.0, 8.0
    mu = U/2  # Half-filled chemical potential.
    beta_list = vcat(2.0 .^ vcat(collect(-15:1:-2), collect(-1.5:0.5:1)),
        collect(3.0:1:16))
    return (;
        application=:SU2SU2_Hubbard_square_lattice_half_filling, symmetry=:su2su2, ensemble=:half_filling,
        L, W, boundary_x=:OBC, boundary_y=W==1 ? :OBC : :PBC, t, U, mu,
        Ne=L*W, Ns=200, beta_list, beta0=1.0,
        settn_bond_dimension=1024, purification_bond_dimension=1024,
        sample_bond_dimension=1024,
        # One-site CBE at every width; D_cbe = D + floor(D*cbe_ratio).
        cbe_ratio=0.25, cbe_tolerance=1e-8, rsvd=true,
        settn_maxorder=6, settn_maxiter=8,
        settn_noise=[(0.1,0.1),(0.1,0.01),(0.1,0.001)],
        settn_tolerance=1e-13, settn_compression=1e-16,
        settn_cutoff=1e-15, purification_cutoff=1e-15, born_cutoff=1e-15,
        sample_cutoff=1e-15, krylov_dimension=32, tdvp_tolerance=1e-12,
        probability_tolerance=1e-10, number_tolerance_per_site=1e-4,
        number_flow_maximum_iterations=9, number_fit_tolerance_per_site=1e-8,
        number_fit_maximum_iterations=128, environment_disk=true,
        gc_step=false, gc_sweep=true,
        # Sampling efficiency eta from Eq. (21); resample below the threshold.
        resample=true, eta_threshold=0.5, eta_weight=:estimator,
        resampling_method=:sir, sir_draw_mode=:with_replacement_until_unique,
        alpha_selection=:loo_multiplicity, resampling_alpha=0.5,
        alpha_max=0.999, alpha_tolerance=1e-8,
        resampling_policy=:threshold_each, forced_betas=Float64[],
        maximum_resampling_events=typemax(Int), eta_target=0.9,
        acceptance_policy=:always,
        thermal_seed=UInt64(20250301), master_seed=UInt64(20250302),
        # Optional checks validate the MPI size and Julia threads set by the launcher.
        expected_mpi_size=nothing, expected_julia_threads=nothing, blas_threads=1,
        output_dir=joinpath(@__DIR__, "output"), reuse_beta0=true, overwrite_output=false,
        save_interval_seconds=3600.0, stop_after_beta=nothing,
        purification=TanTRG(), evolution=TDVP(),
        task=:all,  # :all, :thermal, :purification_reference, :cooling
        measure_local_observables=false,
        dynamics_betas=Float64[],
        dynamics_dt=0.05, dynamics_tmax=20.0,
        frequency_list=collect(0.0:0.05:20.0), window_cutoff=20.0, broadening=0.0,
        ground_state_twice_spin=0, ground_state_dynamics=false,
        ground_state_bond_dimensions=[256,512,1024], ground_state_sweeps=12,
        ground_state_energy_tolerance=1e-9,
    )
end
