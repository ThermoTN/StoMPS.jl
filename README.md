# StoMPS.jl

## Introduction

StoMPS.jl is a package for finite-temperature quamtum many-body calculation, which 
is based on perfect Born sampling of symmetric thermal tensor networks 
([Gao et al., 2026](https://arxiv.org/abs/2609.35754)). 


## Methods

The algorithm combines thermal purification at high temperatures with
stochastic sampling at low temperatures.

### Purification

We first construct thermal density operator
$\rho(\beta/2)=e^{-\beta H/2}/\sqrt{Z}$ in matrix product operator (MPO) form.
We use SETTN (series-expansion thermal tensor network; [Chen et al., 2017](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.161104)) to initialize this MPO at high temperature, 
and then use TanTRG (tangent-space tensor renormalization group; [Li et al., 2023](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.226502)) to evolves it in
imaginary time with the time-dependent variational principle (TDVP;
[Haegeman et al., 2011](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.107.070601),
[Haegeman et al., 2016](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.165116)) to reach
the sampling inverse temperature `beta0`.

### Stochastic sampling

At `beta0`, perfect Born sampling 
sequentially draws auxiliary symmetry-sector configurations according to the
conditional Born probabilities of the purified MPO. Projection onto whole
multiplets preserves non-Abelian symmetries in the resulting stochastic MPS
samples. Each sample is subsequently evolved to lower temperatures using one-site
TDVP with controlled bond expansion (CBE;
[Gleis et al., 2023](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.246402),
[Li et al., 2024](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.026401)).



## A Quick start

Run these commands from the repository root with Julia 1.11 or later:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -t2 -e 'using StoMPS; include("examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/config.quickstart.jl"); run_application(config())'
```

The [quick-start configuration](examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/config.quickstart.jl)
runs a half-filled 2×2 canonical Hubbard calculation with open boundaries
in both directions in one process. See
[setup and launch](docs/src/guide/quickstart.md) for its parameters, resource
measurements, output location and MPI commands. Then use
[reading results](docs/src/guide/results.md#energy-and-completion) to check
completion and print energy per site with its sampling error.



## Examples

These examples cover finite-temperature Hubbard calculations at half filling
and finite doping. The calculations of half-filled Hubbard model on square lattice uses $SU(2) \times SU(2)$ symmetry;
the others use $U(1) \times SU(2)$.

- [Half-filled square-lattice Hubbard model in grand-canonical ensemble](examples/SU2SU2_Hubbard_square_lattice_half_filling/): Thermal observables (e.g. energies, charge density, double occupancy, spin/charge correlation functions) of half-filled square-lattice Hubbard model ($\mu=U/2$) in grand-canonical ensemble.
- [square-lattice Hubbard model with fixed particle number in canonical ensemble](examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/): Thermal observables (e.g. energies, charge density, double occupancy, spin/charge correlation functions) of square-lattice Hubbard model with fixed particle number in canonical ensemble.
- [square-lattice Hubbard model with fixed mean particle number in grand-canonical ensemble](examples/U1SU2_Hubbard_square_lattice_grand_canonical_ensemble/): Thermal observables (e.g. energies, charge density, double occupancy, spin/charge correlation functions) of  square-lattice Hubbard model with fixed mean particle number in grand-canonical ensemble.
- [triangular-lattice Hubbard model with fixed particle number in canonical ensemble](examples/U1SU2_Hubbard_triangular_lattice_canonical_ensemble/): Thermal observables (e.g. energies, charge density, double occupancy, spin/charge correlation functions) of triangular-lattice Hubbard model with fixed particle number in canonical ensemble.
- [Finite-*T* dynamics of 1D Hubbard chain in canonical ensemble](examples/U1SU2_Hubbard_chain_dynamics/): Optical conductivity of 1D Hubbard chain with fixed particle number in canonical ensemble.

See the [application guide](docs/src/guide/applications.md) for configurations
and references.

## Documentation and development

The [manual](docs/src/index.md) links the user guides, method descriptions and
API reference. Follow [convergence and uncertainty](docs/src/guide/validation.md)
when selecting production parameters.

FiniteMPS owns tensor containers and TDVP; StoMPS owns the sampling and shared
driver. Run the Hubbard and triangular Hubbard tests with
`julia --project=. -t2 test/runtests.jl`; see the
[validation guide](docs/src/guide/validation.md#checks-already-available) for its scope.

## Maintainer

- Jianxin Gao, <jxgao@buaa.edu.cn>, Beihang University



## Reference

- J. Gao et al., *Perfect Born Sampling of Symmetric Thermal Tensor Network for Quantum Lattice Models*, [arXiv:2609.35754](https://arxiv.org/abs/2609.35754).
- B.-B. Chen et al., *Series-expansion thermal tensor network approach for quantum lattice models*, [Phys. Rev. B **95**, 161104(R) (2017)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.161104).
- Q. Li et al., *Tangent Space Approach for Thermal Tensor Network Simulations of the 2D Hubbard Model*, [Phys. Rev. Lett. **130**, 226502 (2023)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.226502).
- Q. Li et al., *Thermal tensor network simulations of lattice fermions with fixed filling*, [Phys. Rev. B **113**, 085150 (2026)](https://journals.aps.org/prb/abstract/10.1103/xl6n-9tbh).
- J. Haegeman et al., *Time-Dependent Variational Principle for Quantum Lattices*, [Phys. Rev. Lett. **107**, 070601 (2011)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.107.070601).
- J. Haegeman et al., *Unifying time evolution and optimization with matrix product states*, [Phys. Rev. B **94**, 165116 (2016)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.165116).
- A. Gleis et al., *Controlled Bond Expansion for Density Matrix Renormalization Group Ground State Search at Single-Site Costs*, [Phys. Rev. Lett. **130**, 246402 (2023)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.246402).
- J.-W. Li et al., *Time-Dependent Variational Principle with Controlled Bond Expansion for Matrix Product States*, [Phys. Rev. Lett. **133**, 026401 (2024)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.026401).

