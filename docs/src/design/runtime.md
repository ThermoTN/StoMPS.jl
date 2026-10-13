# Execution, output and extension interfaces

The driver is shared by all examples. `main.jl` initializes MPI and calls
`run_application(config; comm=MPI.COMM_WORLD)`. Library calls may omit `comm`
for serial execution. MPI is an optional extension loaded through `StoMPSMPIExt`.
Rank zero owns the main log and numerical result; ranks own their local MPS
and sample logs. Communication carries scalar records and configurations only.

## Stages

`run_purification(c; comm)` prepares and returns a `ThermalCheckpoint` at beta0.
An optional `before_step(index)` callback can reproduce a reference calculation's
random stream. Record its convention in `c.thermal_seed_schedule`; it participates
in the checkpoint identity. Without a callback, operation-derived seeding is
unchanged. `settn_gc_step` may override `gc_step` for SETTN alone.
`settn_verbose` and `evolution_verbose` optionally enable FiniteMPS sweep output
(both default to zero). SU2SU2 tanTRG retains its environment between endpoints,
as does coupled U1SU2 flow.
With a shared output directory, each MPI rank loads the published file itself.
`run_stochastic(c; checkpoint, checkpoint_path, comm)` samples and cools from
that point. Supply either an in-memory checkpoint on every participating rank
or its shared path. Only rank zero returns the numerical result dictionary.
`output_dir=nothing` supports small serial calculations without file output.

`run_application` dispatches the tasks `:thermal`, `:cooling`, `:all`,
`:purification_reference`, `:ground_state`, `:dynamics` and `:ed`. A full
purification reference stores its numerical curve and the beta0 MPO.
`run_ground_state` saves a converged final MPS; `run_dynamics` loads
only explicitly requested source temperatures. Finite-temperature source
snapshots are recorded after resampling and measurement at that beta.
Independent dynamics runs may select a subset of saved temperatures and change
the real-time grid or Fourier-analysis controls; the thermal source identity
and its scalar records still have to match.

`reuse_beta0=true` checks physical representation, lattice, thermal controls,
thermal seed and preparation grid. Stochastic sample counts, seeds and
resampling controls may change without recomputing that preparation.
The checkpoint retains its original source digest, dependency versions and
physical conventions; stochastic results also retain that preparation provenance.
An existing numerical run is protected unless `overwrite_output=true`.
An interrupted stochastic calculation restarts from beta0. Use a new output
directory or `overwrite_output=true` for the repeated run.

## Algorithms and tensors

`MPS` and `MPO` are FiniteMPS `DenseMPS` subtypes; operators are `SparseMPO`.
Born MPS retain the physical sector-marker leg, including non-Abelian carrier
dimensions.

`EvolutionMethod` dispatches `evolution_step!(state, H, method, delta, config;
rng, seed, D, cutoff)`. The step leaves the state normalized and returns
`log_norm_increment`, `elapsed_seconds`, `maximum_truncation`,
`maximum_entropy` and `maximum_bond_dimension`. `TDVP()` is the current method;
thermal samples pass `delta=-delta_beta/2`. A driver-supplied `seed` preserves
the established FiniteMPS operation stream. Without one, TDVP derives its
default-RNG seed from the explicit `rng`. The driver adds the returned
squared-norm logarithm once to both estimator base `a` and physical norm `h`.

`PurificationMethod` dispatches `thermal_trajectory(method, model, config;
before_step, beta_grid, on_step)` and the endpoint-based `purification_step!`.
`TanTRG()` uses SETTN initialization and one-site TDVP, including coupled H/N
flow for U1SU2. Additional purification methods can define their initialization
and internal temperature schedule through `thermal_trajectory`. Callbacks
receive normalized state and scalar thermodynamic records at the requested
endpoints. DMRG is a separate two-site reference.
Real-time JP propagation dispatches separately on `EvolutionMethod` because its
two propagated vectors must retain their physical amplitudes and energy shift.

## Persistence

The versioned JLD2 format is `stomps-0.6-v1` and rejects earlier schemas.
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
| `dynamics/sources/beta_INDEX/RUN_ID/sample_ID.jld2` | Requested source MPS and matching scalar record |
| `dynamics/sources/beta_INDEX/manifest.jld2` | Complete source population, frozen weights and run identity |
| `results/dynamics_result.jld2` | Numerical JP trajectories and transforms |
| `results/ground_state.jld2` | One final DMRG MPS, sector and convergence status |
| `results/ground_state_result.jld2` | Numerical sweeps and optional ground-state JP |

Numerical results contain no tensors and cannot restore an intermediate
thermal or stochastic state. Dynamics source MPS support the requested JP
calculations; DMRG persists its final MPS. Intermediate DMRG sweep states and
real-time states are not saved. Older research files remain in their existing
format. Generated outputs and local configurations stay Git ignored.

For configurations from before 0.6, replace `resume` with `reuse_beta0` and a
run from beta0. `save_interval_seconds` replaces `checkpoint_interval_seconds`
and controls numerical-history saves.

## Logs and uncertainty

`log/run.log` starts each stage with all effective parameters. A matching
stochastic stage appends to its completed preparation. Root prints fixed-width
tables with ten-decimal scientific notation, `NaN` for unavailable diagnostics,
and a flush after each record. Born and mixture probabilities are formatted
from log probabilities without exponentiating the full probability.

`log/samples/rank_XXX_sample_XXX.log` records individual representatives.
`sample_index` identifies a logging history independently of mutable pool IDs
and RNG lineage. New selected representatives receive new indices; their
lineages and operation seeds retain the original stochastic meaning.
Sample `log_weight` is not ensemble lnZ; physical H-only norm is separate.
Resampling logs explicitly identify pre-event, replay and post-event phases.

MPS/MPO truncation is the squared `TrunErr` from the second TDVP
half-sweep. `time(s)` excludes
observable measurement and resampling; ensemble time is the maximum summed
local evolution time over ranks. Initial beta0 sample time is zero.
Purification diagnostic rows additionally store `cumulative_evolution_seconds`
and `wall_elapsed_seconds`. The latter starts just before the thermal trajectory
and includes initialization, measurement, diagnostics and preceding output I/O;
it excludes model construction and final output after the endpoint is measured.
The launcher must separately measure whole-process wall time when comparing
end-to-end costs, including Julia startup and final file publication.

Energy errors use the multiplicity/fitted-mean estimator.
lnZ/eta errors use a conditional paired delta method: for shifted per-copy
weights w, `mu1=sum(m*w)/Nv`, `mu2=sum(m*w^2)/Nv`, propagate the joint moments
of w and w², with variance denominator `Nv*(Nv-1)`. Hold inherited scale,
adaptive proposals and the current fitted lambda fixed. Undefined eta errors
are NaN.
