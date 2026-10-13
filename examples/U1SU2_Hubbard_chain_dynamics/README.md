# U1SU2_Hubbard_chain_dynamics

U(1) charge × SU(2) spin; an open chain with exact Ne in every Born sample. Saves finite-temperature total JP correlations and a matched ground-state DMRG/connected-JP reference.

`config.jl` defines parameters and `main.jl` initializes MPI. The shared StoMPS
driver owns thermal preparation, sampling, resampling, logging and persistence.

From the repository root:

```sh
julia --project=examples/U1SU2_Hubbard_chain_dynamics -e 'using Pkg; Pkg.instantiate(); using MPI; MPI.install_mpiexecjl()'
mpiexecjl --project=examples/U1SU2_Hubbard_chain_dynamics -n 4 julia -t2 --project=examples/U1SU2_Hubbard_chain_dynamics examples/U1SU2_Hubbard_chain_dynamics/main.jl
```

Edit `config.jl`, or copy it to ignored `config.local.jl`
and append `--config PATH` to the command. Personal output and machine settings
belong in the local copy. Use `task=:thermal` to prepare through beta0,
`:purification_reference` to purify through the full beta grid, `:cooling`
for the sampled trajectory, or `:all` for the complete workflow. MPI size
comes from the launcher; Julia threads from `-t`.
Keep BLAS at one thread unless a measured workload justifies a change.

The package implements the JP calculation. Set `dynamics_betas`, `dynamics_dt`
and `dynamics_tmax`; `task=:dynamics` uses the saved cooling snapshots.
`frequency_list`, `window_cutoff` and `broadening` control the finite-time
Fourier transform with a Parzen window. Keep `window_cutoff=dynamics_tmax`
when changing the time grid. An empty frequency list saves time-domain data only.
Use `task=:ground_state` for matched DMRG/connected JP, or `task=:ed` for
the small-chain fixed-N spin-block reference.

JP dynamics uses an open chain: `W=1`, `boundary_x=:OBC`, `boundary_y=:OBC`.
The first beta is SETTN, purification ends at `beta0`, Born sampling occurs
there, then MPS cooling follows the remaining `beta_list`. All finite-temperature
TDVP steps use one-site CBE with default RSVD and `D_cbe=D+floor(D/4)`.

The supported resampler is with-replacement multiplicity LOO-SIR. Set
`resample=false` for the matched no-resampling control.

Default resampling is enabled, strict estimator `eta < 0.5`, multiplicity
LOO alpha in `[0,0.999]`, and at most one committed event per saved beta.
The default `:with_replacement_until_unique` draws until Ns distinct candidate
MPS IDs or 3Ns attempts, retaining repeated hits as multiplicities. Smaller
ID pools keep all input weights and multiplicities. Select `:with_replacement`
for the fixed-Ns-draw comparison.
`Ns` is the initial count and representative target, `Nv=sum(m)` the actual
virtual population, `K` the stored representatives,
and `Np` the number of different complete sector configurations.

`thermal_seed` controls the thermal parent; `master_seed` controls Born,
resampling and explicitly seeded randomized CBE. Seeds use immutable
lineage/event/step identities, independently of MPI rank. See the
[runtime contract](../../docs/src/design/runtime.md) for MPI ownership.

`output_dir` should be external or ignored. Rank zero writes `log/run.log`;
rank-owned sample logs are `log/samples/rank_XXX_sample_XXX.log`. The numerical
`results/result.jld2` retains every completed beta and event, saved hourly
(`save_interval_seconds=3600.0`) and at stage completion.

`results/checkpoint.jld2` stores the beta0 MPO. `reuse_beta0=true` loads a
compatible preparation; stochastic runs start from beta0. Use a new output
directory or `overwrite_output=true` for a repeated run.
`stop_after_beta` ends the run after saving its numerical history. Requested
`dynamics_betas` save source MPS separately; DMRG saves its final ground-state MPS.
See the shared [output contract](../../docs/src/design/runtime.md) for formats,
restart behavior and uncertainty definitions.

See [physical conventions](../../docs/src/design/physical_conventions.md),
[thermal weights](../../docs/src/design/finite_temperature.md), and
[resampling](../../docs/src/design/resampling.md). Small shared numerical
regressions live in `test/`.
