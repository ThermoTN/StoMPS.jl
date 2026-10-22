# Parameter reference

An application returns a Julia `NamedTuple` from `config()`. Start with its
complete config, including its physical model and numerical settings;
`effective_config` supplies defaults for optional driver controls. Values below describe the shipped
research configs unless marked as driver defaults. The quick start deliberately
uses smaller values. `validate_config(effective_config(c))` checks shared
constraints; geometry, DMRG and measurement checks also run in their stages.

## Model and geometry

Energies share the unit chosen for the hopping `t`; `k_B=1`, so temperature is
an energy and beta is an inverse energy. Setting t=1 expresses energies in
hopping units. See [physical conventions](../design/physical_conventions.md).

| Parameters | Type / usual value | Meaning and constraints |
|---|---|---|
| `application` | `Symbol` | Application label saved in provenance. |
| `symmetry`, `ensemble` | `:u1su2` / `:su2su2`; `:canonical` / `:grand_canonical` / `:half_filling` | Use the combination in the selected application; SU2SU2 is half-filled grand canonical. |
| `lattice_kind` | `Symbol`, driver default `:square` | `:triangular` is explicit and supports U1SU2. |
| `L`, `W` | Positive integers | Physical coordinate lengths; LW sites. A chain has W=1. |
| `boundary_x`, `boundary_y` | `:OBC` or `:PBC` | Independent axes; W=1 requires transverse OBC. SU2SU2 needs a bipartite graph, including even periodic lengths. |
| `t`, `U` | Real energies | Hopping in H has coefficient -t; the sign matters on triangular lattices. |
| `mu` | Real energy, U/2 for SU2SU2 | SU2SU2 energy includes the chemical-potential term. U1SU2 reports the physical Hubbard energy; see [physical conventions](../design/physical_conventions.md). |
| `Ne` | Integer for canonical, real target for grand canonical | Driver requires 0<Ne<2LW. Canonical samples each have Ne; grand canonical fixes the weighted mean. SU2SU2 uses Ne=LW. |

## Temperature and numerical evolution

| Parameters | Type / usual value | Meaning and constraints |
|---|---|---|
| `beta_list`, `beta0` | Real vector; real, usually 1 | Finite, positive, strictly increasing grid. beta0 must be a grid member and greater than its first point. Samples use the grid from beta0 onward. |
| `purification`, `evolution` | `TanTRG()`, `TDVP()` | Shared method objects; finite-T TDVP is one-site CBE. |
| `settn_bond_dimension`, `purification_bond_dimension`, `sample_bond_dimension` | Positive integers | Separate maximum bond cutoffs; assess convergence for each stage. |
| `cbe_ratio`, `cbe_tolerance`, `rsvd` | Nonnegative real, positive real, Bool; 0.25, 1e-8, true | Expansion `D+floor(D*cbe_ratio)`, CBE selection tolerance and randomized SVD. |
| `settn_maxorder`, `settn_maxiter` | Positive integers; 6, 8 | SETTN initialization controls. |
| `settn_noise` | Vector of real pairs | SETTN noise schedule, usually [(0.1,0.1),(0.1,0.01),(0.1,0.001)]. |
| `settn_tolerance`, `settn_compression` | Positive reals; 1e-13, 1e-16 | SETTN convergence and compression controls. |
| `settn_cutoff`, `purification_cutoff`, `born_cutoff`, `sample_cutoff` | Nonnegative reals; 1e-15 | Stage-specific truncation controls passed to tensor truncation; distinct from reported squared TDVP `TrunErr`. |
| `krylov_dimension`, `tdvp_tolerance` | Positive integer, positive real; 32, 1e-12 | Krylov evolution/DMRG controls. |
| `probability_tolerance` | Positive real; 1e-10 | Probability and normalization checks. |
| `number_tolerance_per_site`, `number_flow_maximum_iterations` | Positive real, positive integer; 1e-4, 9 | U1SU2 purification particle-number control; total tolerance scales with LW. |
| `number_fit_tolerance_per_site`, `number_fit_maximum_iterations` | Positive real, positive integer; driver defaults 1e-8, 128 | Grand-canonical fugacity solve; sampled particle numbers must bracket the target mean. |

## Sampling and resampling

| Parameters | Type / usual value | Meaning and constraints |
|---|---|---|
| `Ns` | Integer >1; 100 or 200 | Initial number of Born samples. |
| `thermal_seed`, `master_seed` | Integers, usually `UInt64` | Thermal preparation and stochastic streams derived from operation and sample identifiers. |
| `resample` | Bool; true | Set false for a matched direct-continuation control. |
| `eta_threshold` | Real in [0,1]; 0.5 | Strict trigger eta < threshold. |

The [resampling guide](../design/resampling.md) describes the empirical and
mixture proposals in App. B. Application configurations supply the algorithm controls.

## Measurements and reference tasks

| Parameters | Type / usual value | Meaning and constraints |
|---|---|---|
| `measure_local_observables` | Bool; driver default false | Save site density and double occupancy with errors. |
| `dynamics_betas` | Real vector; driver default empty | Source inverse temperatures for canonical U1SU2 OBC chains; each belongs to the sample grid. |
| `dynamics_dt`, `dynamics_tmax` | Positive step, nonnegative end time | Time units inverse to energy, hbar=1. tmax must be an integer multiple of dt. |
| `frequency_list` | Real vector | Angular frequencies in energy units; an empty vector disables the finite-T frequency transform. |
| `window_cutoff`, `broadening` | Positive time, nonnegative energy | Parzen cutoff must equal the time-grid endpoint for a transform; exponential damping exp(-broadening*time). Need at least two time points. |
| `ground_state_twice_spin` | Nonnegative integer 2S | DMRG spin sector: parity agrees with Ne, and 2S <= min(Ne,2LW-Ne). Thermal canonical sampling includes all spin sectors. |
| `ground_state_bond_dimensions`, `ground_state_sweeps` | Positive integer vector, positive integer | DMRG bond schedule and sweeps per entry. |
| `ground_state_energy_tolerance` | Positive energy; 1e-9 | Successive per-site sweep-energy convergence at the final D. |
| `ground_state_dynamics` | Bool; driver default false | Append connected ground-state JP to the DMRG task; requires a canonical U1SU2 open chain. |

## Execution and output

| Parameters | Type / driver default | Meaning and constraints |
|---|---|---|
| `task` | Symbol, `:all` | See the [task table](../design/runtime.md); DMRG and ED are separate tasks. |
| `output_dir` | String or `nothing`; `nothing` | External or ignored path. MPI preparation and saved-source dynamics require a shared path. Small serial runs with nothing keep results in memory. |
| `reuse_beta0`, `overwrite_output` | Bool; true, false | Compatible parent reuse and explicit replacement of numerical output. Cooling starts from beta0. |
| `save_interval_seconds` | Positive real; 3600 | Numerical save at the next safe boundary, also at completion/failure. |
| `stop_after_beta` | Real or `nothing`; `nothing` | Stop after the first sampled beta >= threshold; a partial curve is incomplete and suppresses automatic dynamics in `:all`. |
| `environment_disk` | Bool; research configs true | Disk contractions require explicit node-local TMPDIR; quick start uses false. |
| `expected_mpi_size`, `expected_julia_threads` | Positive integer or `nothing`; `nothing` | Optional checks against the process and thread counts allocated by the launcher. |
| `blas_threads` | Positive integer; 1 | BLAS threads per process. |
| `gc_step`, `gc_sweep` | Bool | FiniteMPS garbage-collection controls. Optional `settn_gc_step` overrides SETTN alone. |
| `settn_verbose`, `evolution_verbose` | Nonnegative integers; 0 | Optional FiniteMPS sweep output. |
| `thermal_seed_schedule` | Recorded convention, default `:derived` | Set when using an explicit `run_purification(...; before_step)` callback; part of preparation identity. |

When reducing a system, change L/W, Ne, boundary conditions, bond cutoffs and Ns
together. Recheck measurement temperatures, dynamics betas,
and the DMRG spin sector. Keep at least one preparation point before beta0 and
one cooling point after it for a cooling demonstration. Choose a fresh output
directory if thermal preparation parameters changed.
