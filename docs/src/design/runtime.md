# Execution, output and extension interfaces

The driver is shared by all examples. `main.jl` initializes MPI and calls
`run_application(config; comm=MPI.COMM_WORLD)`. Library calls may omit `comm`
for serial execution. MPI is an optional extension loaded through `StoMPSMPIExt`.
Rank zero owns the main log and numerical result; ranks own their local MPS
and sample logs. Communication carries scalar records and configurations.

## Stages

`run_purification(c; comm)` prepares and returns a `ThermalCheckpoint` at beta0.
An optional `before_step(index)` callback can reproduce a reference calculation's
random stream. Record its convention in `c.thermal_seed_schedule`; it participates
in the checkpoint identity. The default uses operation-derived seeding.
`settn_gc_step` may override `gc_step` for SETTN.
`settn_verbose` and `evolution_verbose` optionally enable FiniteMPS sweep output
(both default to zero). SU2SU2 tanTRG retains its environment between endpoints,
as does coupled U1SU2 flow.
With a shared output directory, each MPI rank loads the published file itself.
`run_stochastic(c; checkpoint, checkpoint_path, comm)` samples and cools from
that point. Supply either an in-memory checkpoint on every participating rank
or its shared path. Rank zero returns the numerical result dictionary.
`output_dir=nothing` keeps small serial calculations in memory.

## Tasks

`run_application` selects one task. Paths below are relative to `output_dir`;
numerical files are under `results/`. The application validates its configured
geometry, ensemble and requested measurements before running.

| Task | Prerequisites and stages | Main output |
|---|---|---|
| `:thermal` | Prepare through beta0, or load a compatible parent with `reuse_beta0=true`. | `checkpoint.jld2`; `result.jld2` when preparing a new parent. |
| `:cooling` | Prepare/reuse beta0, then Born sampling and cooling; save requested dynamics sources. | `result.jld2` and optional `dynamics/sources/`. |
| `:all` | Same as cooling; append configured finite-T dynamics when cooling reaches the final sample beta. | Cooling outputs and, conditionally, `dynamics_result.jld2`. |
| `:purification_reference` | Purification through the full `beta_list`. U1SU2 follows mean Ne in the grand-canonical ensemble. | `purification_reference_result.jld2` and beta0 checkpoint. |
| `:dynamics` | Canonical U1SU2 OBC chain; requested complete source manifests and MPS must already exist with matching source identity. | `dynamics_result.jld2`. |
| `:ground_state` | Canonical U1SU2, fixed Ne and total spin; DMRG, optionally followed by connected JP when `ground_state_dynamics=true`. | `ground_state.jld2` and `ground_state_result.jld2`. |
| `:ed` | Canonical U1SU2 OBC chain, L<=8; fixed-particle-number spin-block ED. Memory use grows with the Hilbert-space dimension. | `ed_result.jld2`. |

DMRG, ED and the full purification reference have dedicated tasks.
The return value of `:all` is the cooling result; dynamics has its own result
file. See [reading results](../guide/results.md) for
completion flags and stage success checks.

Finite-temperature source snapshots are recorded after resampling and measurement
at that beta.
Independent dynamics runs may select a subset of saved temperatures and change
the real-time grid or Fourier-analysis controls; the thermal source identity
and its scalar records still have to match.

`reuse_beta0=true` checks physical representation, lattice, thermal controls,
thermal seed and preparation grid. Stochastic sample counts, seeds and
resampling controls can be chosen for each run using that preparation.
The checkpoint retains its original source digest, dependency versions and
physical conventions; stochastic results also retain that preparation provenance.
An existing numerical run is protected unless `overwrite_output=true`.
An interrupted stochastic calculation restarts from beta0. Use a new output
directory or `overwrite_output=true` for the repeated run.

## Algorithms and tensors

`MPS` and `MPO` are FiniteMPS `DenseMPS` subtypes; operators are `SparseMPO`.
Born MPS retain an auxiliary leg for the selected symmetry multiplet, including
its carrier dimension.

`EvolutionMethod` dispatches `evolution_step!(state, H, method, delta, config;
rng, seed, D, cutoff)`. The step leaves the state normalized and returns
`log_norm_increment`, `elapsed_seconds`, `maximum_truncation`,
`maximum_entropy` and `maximum_bond_dimension`. `TDVP()` is the current method;
thermal samples pass `delta=-delta_beta/2`. A driver-supplied `seed` initializes
the FiniteMPS operation stream. With `seed=nothing`, TDVP draws its
default-RNG seed from the explicit `rng`.

`PurificationMethod` dispatches `thermal_trajectory(method, model, config;
before_step, beta_grid, on_step)` and the endpoint-based `purification_step!`.
`TanTRG()` uses SETTN initialization and one-site TDVP, including coupled energy
and particle-number flow for U1SU2. Additional purification methods can define their initialization
and internal temperature schedule through `thermal_trajectory`. Callbacks
receive normalized state and scalar thermodynamic records at the requested
endpoints. DMRG is a separate two-site reference.
Real-time JP propagation dispatches separately on `EvolutionMethod` because its
two propagated vectors must retain their physical amplitudes and energy shift.

## Persistence

The JLD2 format is `stomps-0.6-v1`.
Numerical results contain `parameters`, `observables`, `diagnostics` and
`metadata`. Ordinary histories live in `observables.curve`, pure thermal
histories in `observables.purification`, and replacement audits in
`diagnostics.events`. Metadata records versions, recursive source digest,
geometry, energy convention, seeds, stage completion and progress.

`save_interval_seconds=3600.0` writes at the next safe calculation boundary
after the interval has elapsed. Stage termination also writes.
The last complete ensemble remains separate from in-progress diagnostics.
Temporary files are atomically replaced. Dynamics manifests are published
after all source files are complete, preserving the previous manifest until then.

| File | Contents |
|---|---|
| `results/checkpoint.jld2` | One normalized beta0 MPO and thermal normalization |
| `results/result.jld2` | Numerical thermal/stochastic histories and metadata |
| `results/purification_reference_result.jld2` | Full purification curve |
| `dynamics/sources/beta_INDEX/RUN_ID/sample_ID.jld2` | Requested source MPS and matching scalar record |
| `dynamics/sources/beta_INDEX/manifest.jld2` | Complete source population, frozen weights and run identity |
| `results/dynamics_result.jld2` | Numerical JP trajectories and transforms |
| `results/ground_state.jld2` | One final DMRG MPS, sector and convergence status |
| `results/ground_state_result.jld2` | Numerical sweeps and optional ground-state JP |
| `results/ed_result.jld2` | Small-chain canonical ED reference and optical analysis |

Numerical result files store scalar histories and metadata. Tensor storage
consists of the beta0 checkpoint, requested dynamics source MPS and final
DMRG MPS listed above. Stochastic calculations restart from the beta0
checkpoint. Generated outputs and local configurations stay Git ignored.

## Logs and uncertainty

`log/run.log` starts each stage with all effective parameters. A matching
stochastic stage appends to its completed preparation. Root prints fixed-width
tables with ten-decimal scientific notation, `NaN` for unavailable diagnostics,
and a flush after each record. Born and mixture probabilities are formatted
directly from their logarithms.

`log/samples/rank_XXX_sample_XXX.log` records individual representatives.
`sample_index` identifies a logging history independently of mutable pool IDs
and RNG lineage. New selected representatives receive new indices; their
lineages and operation seeds retain the original stochastic meaning.
Sample `log_weight` stores an estimator log weight; ensemble lnZ and the
physical-H squared-norm logarithm have separate fields.
Resampling logs explicitly identify pre-event, replay and post-event phases.

MPS/MPO truncation is the squared `TrunErr` from the second TDVP
half-sweep. `time(s)` measures the evolution step; ensemble time is the maximum summed
local evolution time over ranks. Initial beta0 sample time is zero.
Purification diagnostic rows additionally store `cumulative_evolution_seconds`
and `wall_elapsed_seconds`. The latter starts just before the thermal trajectory
and includes initialization, measurement, diagnostics and output I/O up to
the endpoint measurement.
The launcher must separately measure whole-process wall time when comparing
end-to-end costs, including Julia startup and final file publication.

Independent runs measure sampling uncertainty. For comparisons and numerical
error estimates, follow the [convergence guide](../guide/validation.md).
