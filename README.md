# StoMPS.jl

Symmetry-preserving Born sampling and finite-temperature Hubbard calculations
with FiniteMPS. This is a Julia 1.11 package, version `0.6.0-DEV`, under
the MIT license. It is not registered in General.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -t2 test/runtests.jl
```

FiniteMPS and FiniteLattices resolve from General. JLD2 is included for results;
MPI is an optional extension. Each example loads this checkout through its
`Project.toml` and supplies parameters to the shared driver.

| Application under `examples/` | Calculation |
|---|---|
| [SU2SU2_Hubbard_square_lattice_half_filling](examples/SU2SU2_Hubbard_square_lattice_half_filling/) | Half filling, grand canonical, mu=U/2; W=1 includes the chain. Energy includes -mu*N. |
| [U1SU2_Hubbard_chain_dynamics](examples/U1SU2_Hubbard_chain_dynamics/) | Open-chain exact-Ne Born sampling, finite-T JP optical response and ground-state reference. |
| [U1SU2_Hubbard_square_lattice_canonical_ensemble](examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/) | Exact-Ne Born sampling, square-lattice energy, density and double occupancy; DMRG reference. |
| [U1SU2_Hubbard_square_lattice_grand_canonical_ensemble](examples/U1SU2_Hubbard_square_lattice_grand_canonical_ensemble/) | Cross-N Born sampling; internal mu fit fixes mean Ne at every saved beta. |
| [U1SU2_Hubbard_triangular_lattice_canonical_ensemble](examples/U1SU2_Hubbard_triangular_lattice_canonical_ensemble/) | Exact-Ne triangular Hubbard; energy and bulk normalized spin correlations at selected temperatures; DMRG reference. |

Every application is configured in `config.jl` and launched through `main.jl`.
For example:

```sh
julia --project=examples/SU2SU2_Hubbard_square_lattice_half_filling -e 'using Pkg; Pkg.instantiate(); using MPI; MPI.install_mpiexecjl()'
mpiexecjl --project=examples/SU2SU2_Hubbard_square_lattice_half_filling -n 4 julia -t2 --project=examples/SU2SU2_Hubbard_square_lattice_half_filling examples/SU2SU2_Hubbard_square_lattice_half_filling/main.jl
```

Use `--config PATH` for a personal ignored config. Defaults use SETTN at the
first beta, mean-N tanTRG purification through beta0=1, Born sampling there,
then MPS cooling with resampling. SU2SU2 fixes mu=U/2 instead of fitting it.
One-site CBE uses RSVD and a 25% expansion. Resampling defaults to estimator
ESS below 0.5 and multiplicity LOO-SIR with `:with_replacement_until_unique`:
draw until Ns distinct candidate MPS IDs or 3Ns attempts, retaining repeated
draws as multiplicities. Pools with fewer than Ns IDs are retained in full
with their original weights. `:with_replacement` retains the fixed-Ns-draw
comparison. The SU2SU2 application also provides mixture, hybrid and WOR methods.

Logs report `Ns` (initial count and representative target), `Nv` (actual virtual
population), `K` (stored representatives) and `Np`
(different complete sector configurations).
Root writes `log/run.log`; each rank writes
`log/samples/rank_XXX_sample_XXX.log`. Tables follow purification, Born sampling,
stochastic evolution and resampling/replay. Timings exclude measurement and
resampling; lnZ/eta errors condition on inherited scale and fitted fugacity.

`results/checkpoint.jld2` stores the beta0 MPO; `reuse_beta0=true` reuses a
compatible preparation. Stochastic runs start from beta0. For a repeated run,
use a new output directory or `overwrite_output=true`.
`results/result.jld2` stores numerical histories every hour
(`save_interval_seconds=3600.0`) and at stage completion. Requested
`dynamics_betas` save source MPS separately, and DMRG saves one final
ground-state MPS. See the [output contract](docs/src/design/runtime.md) for
file contents, restart behavior and format compatibility.

`src/` is organized into `sampling`, `resample`, `measurement`, `evolution`,
`models`, `driver` and `io`.
FiniteMPS owns the `MPS`/`MPO` containers (both `DenseMPS`), tensor operations
and TDVP. Born states retain their rank-four sector-marker tensors.
The sampling APIs use the prefixes `su2su2_*`, `u1su2_*` and `z2_*`.

See [physical conventions](docs/src/design/physical_conventions.md),
[Born sampling](docs/src/design/born_sampling.md),
[thermal weights](docs/src/design/finite_temperature.md),
[resampling methods](docs/src/design/resampling.md) and
[execution, output and extension interfaces](docs/src/design/runtime.md), plus
[contribution boundaries](AGENTS.md).
