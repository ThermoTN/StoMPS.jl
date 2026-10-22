# U1SU2_Hubbard_chain_dynamics

## Model and Hamiltonian

This example computes finite-temperature dynamics of
1D Hubbard chain in the canonical ensemble with $U(1) \times SU(2)$ charge
and spin symmetry. The Hamiltonian is

```math
H=-t\sum_{i=1}^{L-1}\sum_\sigma(c^\dagger_{i\sigma}c_{i+1,\sigma}+\mathrm{h.c.})
  +U\sum_i n_{i\uparrow}n_{i\downarrow}.
```

Born sampling fixes every sample to `Ne` particles, which is preserved during
imaginary- and real-time evolution. Real-time evolution
measures the total JP correlation $C^>_{JP}(\tau)=\langle J(\tau)P(0)\rangle$,
where $P$ is the centered polarization and $J=+i[H,P]$. A finite-time Fourier
transform gives the optical conductivity per chain length. See
[physical conventions](../../docs/src/design/physical_conventions.md#optical-response-and-ground-states)
and [Born sampling](../../docs/src/design/born_sampling.md).

## Parameters

All settings are defined in [config.jl](config.jl). The supplied configuration
uses the following parameters:

| Parameter | Value | Meaning |
|---|---|---|
| `L`, `W` | 32, 1 | Open chain with 32 sites. |
| `boundary_x`, `boundary_y` | `:OBC`, `:OBC` | Open boundaries on both axes. |
| `t`, `U` | 1, 8 | Hopping and onsite repulsion. |
| `Ne` | 28 | Fixed particle number in every sample. |
| `Ns` | 100 | Number of Born samples. |
| `settn_bond_dimension`, `purification_bond_dimension`, `sample_bond_dimension` | 256, 1024, 1024 | Maximum bond dimensions for SETTN initialization, tanTRG purification and sample evolution. |
| `cbe_ratio`, `rsvd` | 0.25, `true` | One-site TDVP uses CBE to $D+\lfloor D/4\rfloor$ and randomized SVD. |
| `resample`, `eta_threshold` | `true`, 0.5 | Resampling is triggered when the estimator sampling efficiency falls below 0.5. |
| `dynamics_betas` | `[1, 2, 4, 8, 16]` | Inverse temperatures at which sample MPS are saved for real-time dynamics. |
| `dynamics_dt`, `dynamics_tmax` | 0.05, 20 | Real-time step and maximum propagation time, in units of $1/t$ with $\hbar=1$. |
| `frequency_list` | `0.0:0.05:20.0` | Angular-frequency grid for the conductivity, in hopping units. |
| `window_cutoff`, `broadening` | 20, 0 | Parzen-window cutoff and exponential damping rate in the finite-time transform. |

Temperatures are specified through the inverse-temperature grid `beta_list`.
With $k_B=1$ and $t=1$, $T=1/\beta$ in hopping units. The default grid is:

```julia
beta_list = vcat(
    2.0 .^ vcat(collect(-15:1:-2), collect(-1.5:0.5:1)),
    collect(3.0:1:16),
)
beta0 = 1.0
```

The five dynamics temperatures are $T/t=1,1/2,1/4,1/8,1/16$. Every value in
`dynamics_betas` must belong to the sampled grid. The default transform uses a
positive Fourier sign and a $1/L$ spatial normalization; its time window and
damping conventions are defined in the
[optical convention](../../docs/src/design/physical_conventions.md#optical-response-and-ground-states).

Start with the shared [quick start](../../docs/src/guide/quickstart.md) and check
[convergence](../../docs/src/guide/validation.md) before choosing production
sample counts, bond dimensions and real-time controls.

The default `environment_disk=true` requires an explicitly configured
node-local `TMPDIR`; see [setup and launch](../../docs/src/guide/quickstart.md).
Copy the configuration to `config.local.jl` for personal runs.

## Tasks and output

Set `task=:all` to prepare or reuse the thermal parent, perform Born sampling,
cool the samples and run finite-temperature JP dynamics at `dynamics_betas`.
Use `:thermal` for preparation only, or `:cooling` followed by `:dynamics` to
run cooling and real-time evolution separately. `:purification_reference`
performs grand-canonical purification at fixed mean `Ne` over the full
`beta_list`. The shared
[task table](../../docs/src/design/runtime.md#tasks) defines stages and outputs.

The separate `:ground_state` task runs DMRG at fixed particle number and spin;
the defaults `ground_state_twice_spin=0` and `ground_state_dynamics=true`
select the singlet and append connected ground-state JP dynamics. The `:ed`
task provides a canonical ED reference for a reduced chain with `L<=8` and a
compatible `Ne`.

Results are written under `output_dir` (the local `output/` directory by
default): `results/checkpoint.jld2` stores the thermal parent at `beta0`, and
`results/result.jld2` stores sampled observables and numerical diagnostics.
Sample MPS for the requested dynamics temperatures are saved under
`dynamics/sources/`, and `results/dynamics_result.jld2` stores JP correlations
and conductivity. Cooling and dynamics have separate completion flags.
The full purification reference uses `results/purification_reference_result.jld2`.
DMRG saves its final MPS in `results/ground_state.jld2` and its numerical
history and connected JP data in `results/ground_state_result.jld2`;
ED uses `results/ed_result.jld2`. See
[reading results](../../docs/src/guide/results.md#jp-dynamics-and-references)
for loading and interpreting these files.

## Related guides

[Applications](../../docs/src/guide/applications.md),
[parameters](../../docs/src/guide/parameters.md),
[resampling](../../docs/src/design/resampling.md),
[results](../../docs/src/guide/results.md),
[convergence](../../docs/src/guide/validation.md),
[troubleshooting](../../docs/src/guide/troubleshooting.md).
