# Applications

Each directory under `examples/` owns its config, environment and MPI entry point.
They all use the shared StoMPS driver. The research configurations are larger
than the [2×2 quick start](quickstart.md).

| Application README | Ensemble and geometry | Distinctive observables / reference |
|---|---|---|
| [SU2SU2 half filling](https://github.com/ThermoTN/StoMPS.jl/blob/main/examples/SU2SU2_Hubbard_square_lattice_half_filling/README.md) | Bipartite square lattice, half-filled grand canonical; W=1 includes chains | Energy includes the chemical-potential term at mu=U/2; whole charge/spin projectors. |
| [U1SU2 chain dynamics](https://github.com/ThermoTN/StoMPS.jl/blob/main/examples/U1SU2_Hubbard_chain_dynamics/README.md) | Open chain, exact Ne in each sample | Total JP correlation and finite-time conductivity; separate fixed-spin DMRG/connected JP and canonical ED tasks. |
| [U1SU2 canonical square lattice](https://github.com/ThermoTN/StoMPS.jl/blob/main/examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/README.md) | Square lattice, exact Ne | Energy, site density and double occupancy; separate DMRG. Contains config.quickstart.jl. |
| [U1SU2 grand-canonical square lattice](https://github.com/ThermoTN/StoMPS.jl/blob/main/examples/U1SU2_Hubbard_square_lattice_grand_canonical_ensemble/README.md) | Square lattice, grand-canonical samples at fitted mean Ne | Energy, density and double occupancy; number-support checks at every saved beta. |
| [U1SU2 canonical triangular lattice](https://github.com/ThermoTN/StoMPS.jl/blob/main/examples/U1SU2_Hubbard_triangular_lattice_canonical_ensemble/README.md) | Triangular lattice, fixed particle number Ne | Energy; separate DMRG. |

U1SU2 reports the physical Hamiltonian H of Eq. (8). Canonical thermal sampling
includes all total-spin sectors, whereas a DMRG run chooses one. Triangular
positive/negative doping are inequivalent at fixed positive hopping. Consult
[physical conventions](../design/physical_conventions.md) before comparing
energies or observables between applications.

Each application README describes its specific configuration. Follow
[setup and launch](quickstart.md) with the chosen
application directory. Public parameter types and coupled constraints live in
the [parameter reference](parameters.md).

Select stages using the [task table](../design/runtime.md), then check
saved outputs and stage completion with [reading results](results.md).
