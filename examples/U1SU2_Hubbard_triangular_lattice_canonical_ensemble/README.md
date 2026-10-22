# U1SU2_Hubbard_triangular_lattice_canonical_ensemble

## Model and Hamiltonian

This example computes finite-temperature properties of the doped
triangular-lattice Hubbard model in the canonical ensemble with
$U(1) \times SU(2)$ symmetry. The Hamiltonian is

```math
H=-t\sum_{\langle ij\rangle,\sigma}(c^\dagger_{i\sigma}c_{j\sigma}+\mathrm{h.c.})
  +U\sum_i n_{i\uparrow}n_{i\downarrow}.
```

Born sampling conditions on `Ne`, so every sample has exactly that particle
number throughout cooling. See
[physical conventions](../../docs/src/design/physical_conventions.md#geometry)
and [Born sampling](../../docs/src/design/born_sampling.md).

## Parameters

All settings are defined in [config.jl](config.jl). The supplied configuration
uses the following parameters:

| Parameter | Value | Meaning |
|---|---|---|
| `lattice_kind` | `:triangular` | Explicit triangular-lattice geometry. |
| `L`, `W` | 12, 4 | Primitive-coordinate columns and rows, with 48 sites. |
| `boundary_x`, `boundary_y` | `:OBC`, `:PBC` | Open along the first primitive direction and periodic along the second, forming a YC cylinder. |
| `t`, `U` | 1, 12 | Hopping and onsite repulsion in Hamiltonian. |
| `Ne` | 56 | Fixed particle number in every sample, corresponding to electron doping $\delta=1/6$. |
| `Ns` | 200 | Initial number of Born samples. |
| `settn_bond_dimension`, `purification_bond_dimension`, `sample_bond_dimension` | 256, 1024, 1024 | Maximum bond dimensions for SETTN initialization, tanTRG purification and sample evolution. |
| `cbe_ratio`, `rsvd` | 0.25, `true` | One-site TDVP uses controlled bond expansion to $D+\lfloor D/4\rfloor$ and randomized SVD. |
| `resample`, `eta_threshold` | `true`, 0.5 | Resampling is triggered when the estimator sampling efficiency falls below 0.5. |
| `bulk_columns` | `5:8` | Central columns used for bulk measurements, containing 16 sites. |
| `correlation_temperatures` | `[0.5, 0.25, 0.125]` | Temperatures at which bulk spin correlations are measured, in hopping units. |

Temperatures are specified through the inverse-temperature grid `beta_list`.
With $k_B=1$ and $t=1$, $T=1/\beta$ in hopping units. The default grid is:

```julia
beta_list = vcat(
    2.0 .^ vcat(collect(-15:1:-2), collect(-1.5:0.5:1)),
    collect(3.0:1:16),
)
beta0 = 1.0
```

Start with the shared [quick start](../../docs/src/guide/quickstart.md) and check
[convergence](../../docs/src/guide/validation.md) before choosing production
sample counts and bond dimensions.

The default `environment_disk=true` requires an explicitly configured
node-local `TMPDIR`; see [setup and launch](../../docs/src/guide/quickstart.md).
Copy the configuration to `config.local.jl` for personal runs.

## Tasks and output

Set `task=:all` to prepare or reuse the thermal parent, perform Born sampling
and cool the samples, measuring bulk spin correlations at the configured
temperatures. Use `:thermal` for preparation only, `:cooling` for
sampling/cooling, or `:purification_reference` for grand-canonical purification
at fixed mean `Ne` over the full `beta_list`. The separate `:ground_state` task
runs DMRG at fixed particle number and spin; `ground_state_twice_spin=0` selects
the singlet sector by default. The shared
[task table](../../docs/src/design/runtime.md#tasks) defines stages and outputs.

Results are written under `output_dir` (the local `output/` directory by
default): `results/checkpoint.jld2` stores the thermal parent at `beta0`, and
`results/result.jld2` stores sampled energies and numerical diagnostics.
At the requested correlation temperatures, a curve row's `correlations` field
contains the normalized spin correlation `css_tot`, its sampling error and
the raw bulk moments. The full purification reference uses
`results/purification_reference_result.jld2`. DMRG saves its final MPS in
`results/ground_state.jld2` and its numerical history in
`results/ground_state_result.jld2`.

## Related guides

[Applications](../../docs/src/guide/applications.md),
[parameters](../../docs/src/guide/parameters.md),
[resampling](../../docs/src/design/resampling.md),
[results](../../docs/src/guide/results.md),
[convergence](../../docs/src/guide/validation.md),
[troubleshooting](../../docs/src/guide/troubleshooting.md).
