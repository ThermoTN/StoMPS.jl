# SU2SU2_Hubbard_square_lattice_half_filling

## Model and Hamiltonian
This example computes finite-temperature properties of the half-filled
square-lattice Hubbard model with $SU(2) \times SU(2)$ charge and spin symmetry. 
The Hamiltonian is 

```math
H=-t\sum_{\langle ij\rangle,\sigma}(c^\dagger_{i\sigma}c_{j\sigma}+\mathrm{h.c.})
  +U\sum_i n_{i\uparrow}n_{i\downarrow}-\mu\sum_i n_i,
```

where $\mu=U/2$. See [physical conventions](../../docs/src/design/physical_conventions.md)
and [Born sampling](../../docs/src/design/born_sampling.md).

## Parameters

All settings are defined in [config.jl](config.jl). The supplied configuration
uses the following parameters:

| Parameter | Value | Meaning |
|---|---|---|
| `L`, `W` | 8, 4 | Cylinder length and circumference, with 32 sites. |
| `boundary_x`, `boundary_y` | `:OBC`, `:PBC` | Open along the length and periodic around the circumference. |
| `t`, `U`, `mu` | 1, 8, 4 | Hopping, onsite repulsion and chemical potential. |
| `Ne` | `L*W=32` | Mean particle number at half filling in the grand-canonical ensemble. |
| `Ns` | 200 | Number of Born samples. |
| `settn_bond_dimension`, `purification_bond_dimension`, `sample_bond_dimension` | 1024 each | Maximum bond dimensions for SETTN initialization, tanTRG purification and sample evolution. |
| `cbe_ratio`, `rsvd` | 0.25, `true` | One-site TDVP uses controlled bond expansion to $D+\lfloor D/4\rfloor$ and randomized SVD. |
| `resample`, `eta_threshold` | `true`, 0.5 | Resampling is triggered when the estimator sampling efficiency falls below 0.5. |

Temperatures are specified through the inverse-temperature grid `beta_list`.
With $k_B=1$ and $t=1$, $T=1/\beta$ in hopping units. The default grid is:

```julia
beta_list = vcat(
    2.0 .^ vcat(collect(-15:1:-2), collect(-1.5:0.5:1)),
    collect(3.0:1:16),
)
beta0 = 1.0
```

Start with the shared
[quick start](../../docs/src/guide/quickstart.md) and check
[convergence](../../docs/src/guide/validation.md) before choosing production
sample counts and bond dimensions.

The default `environment_disk=true` requires an explicitly configured
node-local `TMPDIR`; see [setup and launch](../../docs/src/guide/quickstart.md).
Copy the configuration to `config.local.jl` for personal runs.

## Tasks and output

Set `task=:all` to prepare or reuse the thermal parent, perform Born sampling
and cool the samples. Use `:thermal` for preparation only, `:cooling` for
sampling/cooling, or `:purification_reference` for purification over the full
`beta_list`. The shared
[task table](../../docs/src/design/runtime.md#tasks) defines stages and outputs.

Results are written under `output_dir` (the local `output/` directory by
default): `results/checkpoint.jld2` stores the thermal parent at `beta0`, and
`results/result.jld2` stores sampled observables and numerical diagnostics.
The full purification reference uses `results/purification_reference_result.jld2`.

## Related guides

[Applications](../../docs/src/guide/applications.md),
[parameters](../../docs/src/guide/parameters.md),
[resampling](../../docs/src/design/resampling.md),
[results](../../docs/src/guide/results.md),
[convergence](../../docs/src/guide/validation.md),
[troubleshooting](../../docs/src/guide/troubleshooting.md).
