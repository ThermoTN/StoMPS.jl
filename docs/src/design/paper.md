# Paper and notation

Jianxin Gao et al.,
[*Perfect Born Sampling of Symmetric Thermal Tensor Network for Quantum Lattice
Models*, arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1),
describes the method implemented by StoMPS.

| Quantity | Paper notation | API |
|---|---|---|
| Number of lattice sites | ``N`` | `L*W` |
| Particle number | ``N_{\rm e}`` | `Ne` |
| Initial number of Born samples | ``N_{\rm s}`` | `Ns` |
| Auxiliary configuration | ``\alpha`` | `configuration` |
| Normalized Born sample | ``\lvert\psi_\alpha\rangle`` | `state` |
| Initial Born probability | ``p_{\beta_0}(\alpha)`` | `born_probability` |
| Sampling efficiency | ``\eta`` | `eta` |
| Mixture proportion, App. B | ``\lambda`` | `resampling_alpha` |
| Local charge and spin | ``(C,S)`` or ``(C_z,S)`` | SU2SU2 or U1SU2 sector labels; U1SU2 stores `q=C_z-1` |
| Conductivity per chain length | ``\sigma(\omega)/L`` | `transform.conductivity` |

Section II and Appendix A define [perfect Born sampling](born_sampling.md).
Appendix B presents [annealed Born resampling](resampling.md), Appendix C
describes [imaginary-time evolution](finite_temperature.md), and Appendix D
gives fixed-particle-number sampling. The Hubbard Hamiltonian and optical
response in Secs. III and V use the [physical conventions](physical_conventions.md).
