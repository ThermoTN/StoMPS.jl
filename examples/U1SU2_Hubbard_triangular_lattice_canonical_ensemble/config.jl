# Demonstration parameters; put production settings in ignored config.local.jl.
using StoMPS: TanTRG, TDVP

function config()
    L, W = 12, 4
    t, U = 1.0, 12.0
    beta_list = vcat(2.0 .^ vcat(collect(-15:1:-2), collect(-1.5:0.5:1)),
        collect(3.0:1:16))
    return (;
        application=:U1SU2_Hubbard_triangular_lattice_canonical_ensemble,
        symmetry=:u1su2, ensemble=:canonical, lattice_kind=:triangular,
        L, W, boundary_x=:OBC, boundary_y=W==1 ? :OBC : :PBC, t, U,
        Ne=56, Ns=200, beta_list, beta0=1.0,
        # Physical columns used for bulk measurements.
        bulk_columns=L%3==0 ? (L÷3+1:2*(L÷3)) : nothing,
        # T uses the same energy units as t, k_B=1. Every 1/T must be in beta_list.
        # To measure T=0.4, explicitly add beta=2.5 to beta_list first.
        correlation_temperatures=[0.5,0.25,0.125],
        settn_bond_dimension=256, purification_bond_dimension=1024,
        sample_bond_dimension=1024,
        cbe_ratio=0.25, cbe_tolerance=1e-8, rsvd=true,
        settn_maxorder=6, settn_maxiter=8,
        settn_noise=[(0.1,0.1),(0.1,0.01),(0.1,0.001)],
        settn_tolerance=1e-13, settn_compression=1e-16,
        settn_cutoff=1e-15, purification_cutoff=1e-15, born_cutoff=1e-15,
        sample_cutoff=1e-15, krylov_dimension=32, tdvp_tolerance=1e-12,
        probability_tolerance=1e-10, number_tolerance_per_site=1e-4,
        number_flow_maximum_iterations=9, environment_disk=true,
        gc_step=false, gc_sweep=true,
        resample=true, eta_threshold=0.5, eta_weight=:estimator,
        resampling_method=:sir, sir_draw_mode=:with_replacement_until_unique,
        alpha_selection=:loo_multiplicity, alpha_max=0.999, alpha_tolerance=1e-8,
        maximum_resampling_events=typemax(Int),
        thermal_seed=UInt64(20250301), master_seed=UInt64(20250302),
        expected_mpi_size=nothing, expected_julia_threads=nothing, blas_threads=1,
        output_dir=joinpath(@__DIR__, "output"), reuse_beta0=true, overwrite_output=false,
        save_interval_seconds=3600.0, stop_after_beta=nothing,
        purification=TanTRG(), evolution=TDVP(),
        task=:all, # :all, :thermal, :purification_reference, :cooling, :ground_state
        ground_state_twice_spin=0,
        ground_state_bond_dimensions=[256,512,1024], ground_state_sweeps=12,
        ground_state_energy_tolerance=1e-9,
    )
end
