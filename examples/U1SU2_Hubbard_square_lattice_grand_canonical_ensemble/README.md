# U1SU2_Hubbard_square_lattice_grand_canonical_ensemble

## Model and Hamiltonian

This example computes finite-temperature properties of the doped square-lattice
Hubbard model in the grand-canonical ensemble with $U(1) \times SU(2)$ symmetry. The Hamiltonian is

```math
H=-t\sum_{\langle ij\rangle,\sigma}(c^\dagger_{i\sigma}c_{j\sigma}+\mathrm{h.c.})
  +U\sum_i n_{i\uparrow}n_{i\downarrow}.
```

Born samples span different particle-number sectors,
and each retains its sampled particle number during cooling. The chemical
potential is fitted after Born sampling and at every saved inverse temperature
to make the weighted ensemble mean equal to `Ne`. Energy, local charge density and double occupancy are
measured. See [physical conventions](../../docs/src/design/physical_conventions.md)
and [Born sampling](../../docs/src/design/born_sampling.md).

## Parameters

All settings are defined in [config.jl](config.jl). The supplied configuration
uses the following parameters:

| Parameter | Value | Meaning |
|---|---|---|
| `L`, `W` | 8, 4 | Cylinder length and circumference, with 32 sites. |
| `boundary_x`, `boundary_y` | `:OBC`, `:PBC` | Open along the length and periodic around the circumference. |
| `t`, `U` | 1, 8 | Hopping and onsite repulsion in Hamiltonian. |
| `Ne` | 28 | Mean particle number. |
| `Ns` | 200 | Number of Born samples. |
| `settn_bond_dimension`, `purification_bond_dimension`, `sample_bond_dimension` | 512, 4096, 4096 | Maximum bond dimensions for SETTN initialization, tanTRG purification and sample evolution. |
| `cbe_ratio`, `rsvd` | 0.25, `true` | One-site TDVP uses controlled bond expansion to $D+\lfloor D/4\rfloor$ and randomized SVD. |
| `resample`, `eta_threshold` | `true`, 0.5 | Resampling is triggered when the estimator sampling efficiency falls below 0.5. |
| `measure_local_observables` | `true` | Save local charge density and double occupancy with sampling errors. |
| `number_fit_tolerance_per_site` | `1e-8` | Per-site tolerance for fitting the ensemble mean particle number. |

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
and cool the samples with the mean particle number fixed by the fitted
chemical potential. Use `:thermal` for preparation only, `:cooling` for
sampling/cooling, or `:purification_reference` for grand-canonical purification
at fixed mean `Ne` over the full `beta_list`. The shared
[task table](../../docs/src/design/runtime.md#tasks) defines stages and outputs.

Results are written under `output_dir` (the local `output/` directory by
default): `results/checkpoint.jld2` stores the thermal parent at `beta0`, and
`results/result.jld2` stores sampled observables, fitted chemical potentials
and numerical diagnostics. The full purification reference uses
`results/purification_reference_result.jld2`.

## Related guides

[Applications](../../docs/src/guide/applications.md),
[parameters](../../docs/src/guide/parameters.md),
[resampling](../../docs/src/design/resampling.md),
[thermal weights](../../docs/src/design/finite_temperature.md),
[results](../../docs/src/guide/results.md),
[convergence](../../docs/src/guide/validation.md),
[troubleshooting](../../docs/src/guide/troubleshooting.md).
