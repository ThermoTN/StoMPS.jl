# Physical and lattice conventions

The Hubbard Hamiltonian, symmetry labels and optical response follow
Secs. II, III and V of [arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1).
See [paper and notation](paper.md) for the API mapping.

## Hubbard models

The U1SU2 applications evolve the physical Hamiltonian

```math
H=-t\sum_{\langle ij\rangle,\sigma}(c^\dagger_{i\sigma}c_{j\sigma}+\mathrm{h.c.})
  +U\sum_i n_{i\uparrow}n_{i\downarrow}.
```

Their reported energy is ``\langle H\rangle``. Purification uses a
separate particle-number MPO to prepare the target mean particle number `Ne`.
The canonical applications condition on `Ne` in the Born step.
The grand-canonical application samples across particle-number sectors, then
fits the chemical potential so that their weighted mean equals `Ne` at each
saved beta. See [symmetric Born sampling](born_sampling.md#Symmetric-Born-sampling).

The SU2SU2 half-filled model evolves and reports

```math
H-\mu\sum_i n_i,\qquad \mu=U/2.
```

The additive energy constant is zero. Config specifies mu explicitly.
Energy rows are per site; sample energies used in the estimators are total. SU2SU2 charge
symmetry fixes the mean particle number to `L*W`, while its projectors select the charge and spin
multiplets ``(C,S)=(1/2,0),(0,1/2)`` of Sec. II B.

## Geometry

`L` and `W` are physical x/y lengths. `boundary_x` and `boundary_y` accept
`:OBC` and `:PBC` independently. `hubbard_model(c; lattice_kind=:square)`
constructs the physical square lattice, then applies the snake ordering.
Each undirected nearest
neighbor is counted once, with both hopping directions included.
Width-two OBC/PBC have the same unique bonds but different declared geometry.
W=1 uses transverse OBC. SU2SU2 requires a bipartite graph; periodic axes
must be even. Sublattice signs derive from physical coordinates.
The public factor order is charge SU(2), then spin SU(2).


## U1SU2 basis

The public local order is empty, single and double occupation, with
`(q,S)=(-1,0),(0,1/2),(1,0)`, charges `q=n-1` and carrier dimensions 1,2,1.
The full spin doublet is retained in the single-occupation sector.
Canonical Born sampling conditions the entire sequential law on `Ne`,
using the remaining-charge projectors in App. D. Grand-canonical sampling
allows different sample particle numbers and fits their weighted mean.
Each sample retains its particle number during evolution; resampling can
change the population of particle-number sectors in the grand-canonical case.

## Optical response and ground states

The chain uses OBC, centered positions, $P=\sum_j x_j n_j$ and $J=+i[H,P]$.
Write real time as $\tau$ to distinguish it from hopping $t$.
It saves the total extensive greater correlation
$C^>_{JP}(\tau)=\langle J(\tau)P(0)\rangle$.
The mixed-bra contraction includes its normalization coefficient and uses
serial current-tree evaluation. The $\tau=0$ result is compared with direct MPO
application. Frequency-transform controls apply to analysis of the saved
time-domain data. Ground-state JP uses $(P-\langle P\rangle)|0\rangle$ and
the $H-E_0$ phase convention.

With a perturbation $H'(\tau)=-F(\tau)P$, linear response gives the causal kernel

$$
K_{JP}(\tau)=i\langle[J(\tau),P]\rangle
           =-2\mathrm{Im}C^>_{JP}(\tau),\qquad \tau\geq0.
$$

The second equality uses Hermiticity of $J$ and $P$. The saved conductivity
describes the response of current **per chain length**, with the following Fourier convention:

$$
\frac{\sigma(\omega)}{L}=\frac1L\int_0^{\tau_{\max}}d\tau\,
 e^{+i\omega\tau} W(\tau/\tau_c)e^{-\gamma\tau}K_{JP}(\tau),
\qquad
W(x)=\begin{cases}
1-6x^2+6x^3,&0\leq x\leq1/2,\\
2(1-x)^3,&1/2<x\leq1,\\
0,&x>1.
\end{cases}
$$

Here $\tau_c$ is `window_cutoff`, $\gamma$ is `broadening` and $\omega$ comes from
`frequency_list`. `OpticalED.finite_time_conductivity` implements this transform
for both ED and sampled dynamics with the composite trapezoidal rule on the
saved time grid. It requires $\tau_c$ to equal the grid endpoint; at least two
time points are needed. Damping uses $\gamma$; sampling efficiency is $\eta$.
The maximum propagation time and window determine the frequency resolution,
as described in Sec. V.

The code uses $k_{\mathrm B}=\hbar=1$, unit charge and unit lattice spacing
$a=1$. Energies and angular frequencies share the hopping energy unit; real
time and $\beta$ use its inverse. Restoring charge and spacing, this
one-dimensional conductivity per chain length has unit $e^2 a/\hbar$.
`Cgreater_JP` is the total correlation, and `transform.conductivity` is
$\sigma(\omega)/L$, matching the ordinate in Fig. 5.

In finite-temperature output, `transform.kernel` is $K_{JP}(\tau)$,
`transform.damping` is $W(\tau/\tau_c)e^{-\gamma\tau}$, and
`transform.conductivity` is complex $\sigma(\omega)/L$, ordered as
`frequency_list`. An empty `frequency_list` leaves `transform=nothing`. Ground-state
output stores connected time-domain JP for subsequent analysis.
See [reading results](../guide/results.md) for the paths and time-domain errors.

Propagation can be checked by comparing canonical MPS JP with Fock-space ED
on the Hubbard dimer. For the transform, a two-level
system with $H=\mathrm{diag}(0,\Delta)$ and off-diagonal elements
$P_{01}=P_{10}=1$ has

$$
K_{JP}(\tau)=2\Delta\tanh(\beta\Delta/2)\cos(\Delta\tau).
$$

Its analytic windowed Laplace transform provides a reference for both real and
imaginary $\sigma$, including the sign, length factor and nonzero broadening.

DMRG uses a specified particle-number and total-spin sector and two-site sweeps. A
reported sweep convergence concerns a chosen $D$ and spin sector. Convergence
in bond dimension and across spin sectors is assessed separately.
Local density and double occupancy are measured
on the same model/geometry as the thermal application.
