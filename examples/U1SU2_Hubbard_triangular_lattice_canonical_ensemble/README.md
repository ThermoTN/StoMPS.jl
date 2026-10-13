# Canonical triangular Hubbard example

This application calculates the doped triangular Hubbard model with charge
U(1) and spin SU(2). Every Born sample has exactly `Ne` particles. It reports
the physical H-only energy per site at every sampled beta and the manuscript
bulk spin correlation at selected temperatures. It includes an optional
fixed-(Ne,S) two-site DMRG reference.

## Configuration and launch

The illustrative defaults in [config.jl](config.jl) are a 12 by 4 YC cylinder
(x open, y periodic), t=1, U=12, Ne=56, Ns=200, and bulk columns 5:8.
The SETTN, purification and sample bond cutoffs are 256, 1024 and 1024.
Choose bond cutoffs and sample counts using convergence checks for the target
calculation.

From the repository root:

```sh
julia --project=examples/U1SU2_Hubbard_triangular_lattice_canonical_ensemble -e 'using Pkg; Pkg.instantiate(); using MPI; MPI.install_mpiexecjl()'
mpiexecjl --project=examples/U1SU2_Hubbard_triangular_lattice_canonical_ensemble -n 4 julia --project=examples/U1SU2_Hubbard_triangular_lattice_canonical_ensemble -t2 examples/U1SU2_Hubbard_triangular_lattice_canonical_ensemble/main.jl
```

For `environment_disk=true`, the launcher must set `TMPDIR` to an existing
node-local scratch directory on every node. Use the MPI backend configured
for this environment; a custom system MPI requires an ABI check and two-rank
smoke before production.

Copy `config.jl` to ignored `config.local.jl`, edit it, and pass
`--config examples/U1SU2_Hubbard_triangular_lattice_canonical_ensemble/config.local.jl`.
Output must be external or Git ignored. `L` and `W` count primitive-coordinate
columns and rows; the vectors are `(sqrt(3)/2, 1/2)` and `(0, 1)`.
Both boundary conditions are configurable. Ne below/above LW gives negative/
positive doping `delta=Ne/(LW)-1`; these are inequivalent for fixed positive t.

`bulk_columns` is an explicit contiguous physical column range. The default
is the central third only when L is divisible by three; otherwise set it
explicitly. Changing L or W also requires choosing an appropriate Ne.

`correlation_temperatures=[0.5,0.25,0.125]` requests beta=2,4,8. T uses the
same energy units as t, with k_B=1; output includes T and T/t. Each inverse
temperature must already occur in `beta_list` at or after `beta0`. For T=0.4,
explicitly add beta=2.5 to the grid. Missing or duplicate points fail validation.
Use `Float64[]` to disable correlation measurements.

## Observable and output

For bulk B and its unique nearest-neighbor bonds with both endpoints in B,

```math
Q_B=\sum_{\{i,j\}\subset B,\ {\rm NN}}\mathbf S_i\cdot\mathbf S_j,
\qquad
\bar p_B=\frac{\langle\sum_{i\in B}n_i\rangle
 -2\langle\sum_{i\in B}n_{i\uparrow}n_{i\downarrow}\rangle}{N_B},
\qquad
\mathcal C_{\rm ss}^{\rm tot}=\frac{4\langle Q_B\rangle}{9N_B\bar p_B^2}.
```

The factor 1/3 follows from SU(2) invariance; the factor 4 comes from the
experimental longitudinal-spin convention. The default bulk has 16 sites
and 40 retained bonds, but the spatial normalizer is **48**, not 40. The
denominator is the squared bulk site-averaged singles probability, not
`<p_i*p_j>` and not a sample-wise normalization. The derivation and comparison
with Lebrat et al., Nature 629, 317 (2024), are in
[physical conventions](../../docs/src/design/physical_conventions.md).

`results/result.jld2` retains all beta rows with energy, weights, particle counts
Ns/Nv/K/Np, configurations, norm diagnostics and resampling events. At requested
temperatures, `row.correlations` contains `css_tot`, `css_tot_stderr`,
`spin_sum`, `bulk_single_mean`, `bulk_density`, `bulk_doping`, their moment
errors and normalization metadata. At other beta points it is `nothing`.
`row.correlations.raw` saves the three scalars `spin_sum`, `number` and `double`
per stored representative, keyed by `id` and aligned with `row.scalars` and
`row.log_weight`. No per-bond measurements are saved.
In `row.scalars`, `h` is the accumulated physical H-only log norm and `a`
is the estimator log weight. Sampled SIR branches reset `a` while retaining `h`.

The default `sir_draw_mode=:with_replacement_until_unique` draws iid candidate
IDs until Ns distinct MPS IDs or 3Ns attempts. Repeats retain multiplicities;
smaller pools keep all original weights and multiplicities. Ns is the initial
count and MPS target, Nv=sum(m) is the actual virtual population, K counts stored
MPS, and Np counts distinct full configurations. Statistical normalization uses
Nv. The `:with_replacement` option retains fixed-Ns-draw SIR.

The ratio is formed after multiplicity-weighted ensemble averaging. Its
first-order error includes numerator/denominator covariance. Bulk density can
differ from the globally fixed density. Measurements occur after SIR at the
requested beta.

`results/result.jld2` records numerical history, source digest, versions,
configuration, geometry, energy convention, sector order and measurement mapping.
It is atomically saved hourly and at completion. A failed measurement preserves
the last complete curve and records its raw diagnostics.

`results/checkpoint.jld2` stores the beta0 MPO. `reuse_beta0=true` reuses
compatible preparation; stochastic runs start at beta0.
`stop_after_beta` ends with numerical results; use a new output directory or
explicit `overwrite_output=true` to repeat from beta0.
DMRG stores one final MPS. Root writes `log/run.log`, ranks write
`log/samples/rank_XXX_sample_XXX.log`. See the shared
[output contract](../../docs/src/design/runtime.md). Seeds remain independent of rank.

## Optional tasks and validation

`task=:all` prepares the beta0 parent and runs canonical cooling. `:thermal`
prepares only the parent; `:cooling` loads it or prepares it if absent.
`:purification_reference` evolves the mean-Ne grand-canonical purification
through the full grid; this is a separate ensemble reference, not a canonical
spin-correlation run. `:ground_state` runs only DMRG with
`ground_state_twice_spin=2S`, whose parity must agree with Ne. DMRG uses the
same Hamiltonian and geometry; sweep convergence alone does not establish
convergence in bond dimension or across total-spin sectors.

Core tests independently construct fixed-N Fock-space Hamiltonians with
fermion signs on a four-site cluster containing two triangles. They compare
energy, physical spin, singles and normalized correlations for Ne=3 and 5.
Shared workflow checks cover covariance propagation, consecutive SIR events,
MPI n=1/2, beta0 reuse and a DMRG energy reference.
