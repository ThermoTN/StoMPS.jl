# Perfect Born sampling

StoMPS samples normalized purification MPOs sequentially in a symmetry-sector
basis, following Sec. II and App. A of
[arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1).
The [notation guide](paper.md) relates the API to the paper. The shared projector loop lives in
`src/sampling/sector_born.jl`; backend files own their distinct local
space, retained auxiliary multiplets, and result fields.
The supported backends are SU(2)×SU(2) and U(1)×SU(2).

## Symmetric Born sampling

At `beta0`, Born sampling sequentially selects whole irreducible multiplets
of the auxiliary local space. Their conditional probabilities are obtained
from the squared norms of the projected tensors. The selected multiplet is
retained as an auxiliary leg of the normalized sample, preserving its full
carrier space. The exact projection uses `notrunc()`.

### SU(2)×SU(2)

The symmetry basis in `su2su2_born_space()` is labeled by charge pseudospin
and spin, in that order. The four local Hubbard states form two doublets:

| Multiplet | Charge and spin `(C,S)` | Local states | Multiplet dimension |
|---|---|---|---|
| Charge | `(1/2,0)` | Empty and doubly occupied | 2 |
| Spin | `(0,1/2)` | Singly occupied, spin up and spin down | 2 |

The projectors select these complete doublets, as in
Sec. II B. Sampling resolves the multiplet label while retaining both states
within the selected doublet, preserving charge and spin SU(2) symmetry.

### U(1)×SU(2)

The symmetry basis in `u1su2_born_space()` is labeled by U(1) charge `q=n-1`
and spin `S`. The public sector order is empty, single and double occupation:

| Multiplet | Charge and spin `(q,S)` | Local states | Multiplet dimension |
|---|---|---|---|
| Empty | `(-1,0)` | Empty | 1 |
| Single | `(0,1/2)` | Singly occupied, spin up and spin down | 2 |
| Double | `(1,0)` | Doubly occupied | 1 |

The single-occupation projector retains the full spin doublet. Each sampled
configuration has a definite total particle number, given by the sum of its
local occupations. This basis supports two ensemble treatments.

#### Canonical ensemble: fixed-particle-number sampling

`u1su2_fixed_n_purification_born_sample` conditions the entire Born rule on
the target particle number `Ne`. At each site, the conditional probability
includes only suffix configurations whose particle number completes that
target, following App. D, Eqs. (24)–(25). Every resulting sample has exactly
`Ne` particles and remains in that sector during imaginary-time evolution,
giving a canonical-ensemble calculation.

#### Grand-canonical ensemble: fitting the mean particle number

`u1su2_purification_born_sample` draws configurations from the full Born
distribution without conditioning on a target particle number. Individual
samples retain their sampled particle numbers, while the ensemble spans
different particle-number sectors.

After sampling, the chemical potential is fitted so that the weighted mean
particle number equals the target `Ne`. The fit is repeated at every saved inverse temperature, including
`beta0`, and after resampling. Thus the grand-canonical calculation fixes the
ensemble mean while allowing particle-number fluctuations. See
[thermal evolution and weights](finite_temperature.md#Grand-canonical-weights)
for the weight convention.

Both treatments start from a U(1)×SU(2) purification prepared at the target
mean `Ne`. Fixed-number conditioning or fitting the sampled mean then
determines the ensemble used for the stochastic calculation.

## Public API

Born and resampling draws use an explicit `rng::AbstractRNG`.
`continue_born_sample!` passes that RNG to its `step!` callback, which owns
any default-RNG seeding needed by FiniteMPS. The shared driver derives seeds
from thermal or master seeds and operation/lineage identities. For TDVP,
`evolution_step!` seeds Julia's default RNG from its `seed` keyword, or draws
a seed from `rng` when that keyword is omitted. See [runtime.md](runtime.md).

`continue_born_sample!` returns named tuples with exactly `beta`, `observable`,
and `relative_log_weight`. It owns beta validation, the evolution-step and
normalization policy; callers provide `step!` and `observe`.

`logweighted_mean(values, log_weights)` is the stable ratio estimator.
`estimate_logpartition(logZ0, relative_log_weights)` adds `logZ0` exactly once
to the internal log-mean-exp. Inputs are nonempty, length-compatible, finite,
and the internal `_logmeanexp` remains private.
