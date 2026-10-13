# Born sampling contract

StoMPS samples normalized purification MPOs sequentially in a symmetry-sector
basis. The shared projector loop lives in
`src/sampling/sector_born.jl`; backend files own their distinct local
space, marker, and result contracts.

## Measurements

The Z₂ backend uses rank-one local sector markers. The public constructors are
`z2_born_space(0,1)` for the physical space and `z2_born_space(1)` for the odd
auxiliary space.

The ``SU(2)\times SU(2)`` backend uses the degenerate projectors `{P_eta, P_spin}` on exactly
`su2su2_born_space()`. Both admitted irreps have multiplicity one but carrier
dimension two. The sampler retains the selected carrier as the third rank-4
marker leg and never splits it into magnetic quantum numbers.

The fixed-N Hubbard backend uses exactly `u1su2_born_space()` and the ordered
projectors `{P_empty,P_single,P_double}`. Its U(1) charge is `q=n-1`, so the
ordered sectors are `(-1,0)`, `(0,1/2)`, `(1,0)` and the complete carrier
dimensions are `[1,2,1]`. It measures the full spin doublet and never resolves
`m_S`.

The parent purification may be grand canonical, but the sampling law is
conditioned on one target particle number:

```text
p_beta0,N(X) = p_GCE(X) / Pr_GCE(N),   sum_site n_site(X) = N.
```

The implementation contracts suffix transfer environments grouped by the
remaining U(1) charge. The selected local conditional includes only suffixes
that can reach `target_N`; this is not rejection sampling. Public result fields
separate `born_probability=p_beta0,N(X)`,
`unconditioned_born_probability=p_GCE(X)`, and
`target_sector_weight=Pr_GCE(N)`.

At each site, the squared norm of each projected candidate is its conditional
Born probability. Candidates sum to one up to the reported numerical residual.
After selection, normalize the projected tensor before moving the canonical
center. The exact snapshot convention is `notrunc()`.

## Unconditioned U(1)×SU(2) law

`u1su2_purification_born_sample` and `u1su2_purification_born_project` use the
same ordered complete projectors, without the fixed-N suffix constraint.
They reuse the shared sequential projector loop. For `P_X` the tensor product
of the selected marker embeddings, `p0(X)=norm(A0*P_X)^2`; completeness of the
empty/single/double projectors gives `sum_X p0(X)=1`. The spin doublet is one
rank-two projector, so its carrier multiplicity is included in the norm.
No additional spin degeneracy factor belongs in an estimator.

Each resulting physical state has `particle_number=sum(occupancies)` and
`total_charge=particle_number-L`, even though the overall sample population
spans N. The rank-4 marker legs retain the selected whole carriers; these are
purified snapshots and do not impose a physical total-spin singlet. The
unconditioned `log_probability` and `born_probability` refer to `p0(X)`.
They obey `p0(X)=Pr(N_X)*p0,N_X(X)`.

## Result diagnostics

A Born result retains:

- the normalized snapshot;
- sampled configuration/sector indices;
- the exact sector ordering used for interpretation;
- local conditional probabilities;
- their summed log probability;
- post-measurement compression norm; and
- maximum probability-normalization residual.

``SU(2)\times SU(2)`` and U(1)×SU(2) results additionally retain `marker_dimensions` in their
documented sector order. Fixed-N results also retain occupancies, target N,
target charge, and conditional and unconditioned probabilities. Changes to
these fields are public result-schema changes.

## Sampling versus enumeration

Production Monte Carlo samples are already distributed according to the Born
probability. Their thermal continuation estimator uses only relative norm
weights; do not multiply by `born_log_probability` again.

The deterministic ``SU(2)\times SU(2)`` four-site oracle visits all 16 configurations
uniformly, so it must include `log(p_configuration)` when combining those
trajectories. That correction belongs only to the private uniform-enumeration
oracle, not to production Monte Carlo.

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

## Required regressions

Z₂ tests enumerate a complete small input basis, check identity and probability
normalization, compare sequential and direct-column probabilities, and compare
each collapsed snapshot with a private direct contraction. ``SU(2)\times SU(2)`` tests check
projector identity/orthogonality, direct/sequential equality, retained carrier
markers, deterministic repeatability, and the two-site ED energy oracle.
Fixed-N U(1)×SU(2) tests enumerate every two-site occupancy configuration in
all five particle-number sectors, compare sequential and direct conditional
probabilities, check the `[1,2,1]` markers and fixed total charge, and compare
each resampling local RDM with a full-chain projector contraction.
Unconditioned U(1)×SU(2) tests enumerate the complete two-site configuration
space, check `p0(X)=Pr(N_X)*p0,N_X(X)`, the identity sector distribution
`[1,4,6,4,1]/16`, marker spaces, normalization and seeded repeatability.
The grand-canonical application additionally compares complete four-site
Fock-space thermodynamics and all 81 sector configurations with ED.
