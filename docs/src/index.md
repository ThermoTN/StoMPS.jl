# StoMPS.jl

```@meta
CurrentModule = StoMPS
```

Symmetry-preserving sequential Born sampling, Hubbard thermal/sample kernels,
stable estimators, and a shared serial/MPI calculation driver.

## Maintained applications

- `examples/SU2SU2_Hubbard_square_lattice_half_filling`: grand-canonical
  half filling, mu=U/2, energy includes -mu*N; W=1 is the chain case.
- `examples/U1SU2_Hubbard_chain_dynamics`: exact-Ne open-chain finite-T JP
  response and its matched ground-state DMRG/connected-JP reference.
- `examples/U1SU2_Hubbard_square_lattice_canonical_ensemble`: exact-Ne Born
  samples, energy and local observables, fixed-sector DMRG reference.
- `examples/U1SU2_Hubbard_square_lattice_grand_canonical_ensemble`: cross-N
  Born samples with internal mu fitting to mean Ne at every saved beta.
- `examples/U1SU2_Hubbard_triangular_lattice_canonical_ensemble`: exact-Ne
  triangular energy, bulk spin correlations at selected T, and DMRG reference.

Every application has its own `config.jl` and `main.jl`; MPI, output and
numerical workflows are shared through StoMPS.
Logs distinguish target Ns, actual virtual count Nv, stored
representatives K and distinct sector configurations Np.

## Numerical contracts

- [Physical and lattice conventions](design/physical_conventions.md)
- [Tensor and symmetry conventions](design/tensor_conventions.md)
- [Born sampling](design/born_sampling.md)
- [Thermal evolution and weights](design/finite_temperature.md)
- [Execution, output and method interfaces](design/runtime.md)
- [Resampling methods](design/resampling.md)
- [Purification/Born derivation](finite_temperature_resampling.md)

## Library API

```@autodocs
Modules = [StoMPS, StoMPS.SU2SU2Fermion]
```
