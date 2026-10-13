# Physical and lattice conventions

## Hubbard models

The U1SU2 applications evolve the physical Hamiltonian

```math
H=-t\sum_{\langle ij\rangle,\sigma}(c^\dagger_{i\sigma}c_{j\sigma}+\mathrm{h.c.})
  +U\sum_i n_{i\uparrow}n_{i\downarrow}.
```

Their reported energy excludes chemical potential. Purification uses a
separate number MPO to follow a grand-canonical fugacity at mean Ne. The
canonical applications condition exactly on Ne only in the Born step.
The grand-canonical application samples across N and fits lambda=beta*mu.

The SU2SU2 half-filled model instead evolves and reports

```math
K=H-\mu N,\qquad \mu=U/2.
```

There is no additive U/4 per-site shift. Config specifies mu explicitly.
Energy rows are per site; sample energies used in the estimators are total. SU2SU2 charge
symmetry fixes mean N=L*W, while its projectors measure eta/spin sectors,
not definite local occupations or magnetic sublevels.

## Geometry

`L` and `W` are physical x/y lengths. `boundary_x` and `boundary_y` accept
`:OBC` and `:PBC` independently. `hubbard_model(c; lattice_kind=:square)`
constructs the physical square lattice, then applies the snake ordering.
Each undirected nearest
neighbor is counted once, with both hopping directions included.
Width-two OBC/PBC have the same unique bonds but different declared geometry.
W=1 uses transverse OBC. SU2SU2 requires a bipartite graph; periodic axes
must be even. Sublattice signs derive from physical coordinates, not the
snake index. The public factor order is charge SU(2), then spin SU(2).

The explicit option `lattice_kind=:triangular` supports U1SU2 and uses unit
primitive vectors `(sqrt(3)/2, 1/2)` and `(0, 1)`. Here x/y label integer
primitive coordinates, not Cartesian components. Positive bond directions
can be chosen as `(1,0)`, `(0,1)` and `(1,-1)`. YC means open x and periodic y.
The model has no A/B sublattice. The square default, hopping sign `-t` and
energy convention are unchanged. On a triangular lattice, positive and
negative doping are physically inequivalent at fixed positive t.

## U1SU2 basis

The public local order is empty, single and double occupation, charges
q=n-1 and carrier dimensions 1,2,1. Fixed-N Born sampling conditions the
entire sequential law on Ne; it does not reject afterward or adjust sample
occupations. Every resulting state and every continued/resampled state
remains in that exact N sector. Unconditioned Born draws retain their N.

## Bulk normalized spin correlation

The triangular canonical application implements the manuscript's spatial
restriction of the normalized spin correlator in Lebrat et al.,
[Nature 629, 317 (2024)](https://www.nature.com/articles/s41586-024-07272-9),
equations (2), (5) and (9) in the
[full text](https://arxiv.org/html/2308.12269). Its normalization is the square
of a spatially averaged singles probability, not a joint singles probability.

The local spin and singles projector follow directly from fermionic number
operators:

```math
S_i^z=\frac{n_{i\uparrow}-n_{i\downarrow}}{2},\qquad
p_i=(n_{i\uparrow}-n_{i\downarrow})^2
    =n_i-2n_{i\uparrow}n_{i\downarrow}.
```

Empty and doubly occupied states have zero spin. The fixed-N thermal trace
includes all total-spin sectors and their magnetic degeneracies; it does not
fix total spin or magnetization. For this SU(2)-invariant ensemble,
`<S_i^a>=0` and `<S_i^a S_j^b>=delta_ab*<S_i dot S_j>/3`. Therefore the
connected longitudinal correlator is exactly one third of the scalar spin
correlator. Merely using an SU2 tensor representation would not justify this
replacement for a magnetized or magnetic-sublevel-selected ensemble.

For bulk sites B and unique undirected nearest-neighbor bonds E_B with both
endpoints in B, define

```math
N_B=|B|,\quad Q_B=\sum_{\{i,j\}\in E_B}\mathbf S_i\cdot\mathbf S_j,
\quad \bar p_B=\frac{1}{N_B}\sum_{i\in B}\langle p_i\rangle.
```

Then

```math
\mathcal C_{\rm ss}(i,j)
 =\frac{4\langle S_i^z S_j^z\rangle_c}{\bar p_B^2}
 =\frac{4\langle\mathbf S_i\cdot\mathbf S_j\rangle}{3\bar p_B^2},
\qquad
\mathcal C_{\rm ss}^{\rm tot}
 =\frac{4\langle Q_B\rangle}{9N_B\bar p_B^2}.
```

`hubbard_spin_correlation_mpo(nsites, bonds)` builds Q_B with physical spin
normalization. Its input must contain distinct undirected, non-self bonds;
it includes no factor 4/3 and no spatial normalization. The two-spin spectrum
is -3/4 in the singlet and +1/4 in the triplet, with zero on spinless states.
These identities and a separate fixed-N Fock-space ED validate the MPO and
the thermal observable for both electron and hole doping.

Bulk columns are selected by physical primitive coordinates before filtering
bonds, independently of snake ordering. The default 12 by 4 YC geometry with
columns 5:8 has N_B=16 and 40 internal bonds (16,12,12 by direction). The
normalizer remains 3N_B=48. Dividing by 40 would instead define a retained-bond
average and change the requested observable. If all retained bond correlators
were identical, the specified total would be 5/6 of that correlator.

The experiment selects anchor sites by local density and sums three neighbor
directions; the second endpoint need not belong to that selection. This
application selects a fixed physical bulk and requires both endpoints inside.
Its denominator uses the same bulk site average. It reports bulk density
and doping separately from global Ne/(LW)-1. Finite canonical cylinders and
experimental trapped, density-binned snapshots need not yield the same value
even after spatial selection is matched. The spin observable retains its
meaning; the ensemble, spatial weighting and density distribution differ.
No conditional-spin normalization, per-bond density normalization or imaging
correction is applied. A vanishing bulk singles denominator is an error.

## Optical response and ground states

The chain uses OBC only, centered positions, P=sum_j x_j*n_j and J=+i[H,P].
It saves the total extensive greater correlation C_JP(t)=<J(t)P(0)>.
The existing mixed-bra coefficient correction and serial current-tree
contraction are retained. The t=0 result is compared with direct MPO
application. Frequency-transform controls apply to analysis of the saved
time-domain data. Ground-state JP uses (P-<P>)|0> and the H-E0 phase convention.

DMRG uses a specified fixed-N total-spin sector and two-site sweeps. A
reported sweep convergence at a chosen D does not establish convergence in
D or across spin sectors. Local density and double occupancy are measured
on the same model/geometry as the thermal application.
