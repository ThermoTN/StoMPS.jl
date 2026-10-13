# StoMPS.jl Project Instructions

This file governs the repository. Stable numerical contracts live in
[`docs/src/design/`](docs/src/design/); application parameters live beside each
application. Historical profile/server notes are local material, not contracts.

## Workflow

- Read this file before edits. Preserve unrelated user work; discover with `rg`
  and edit manually with `apply_patch`. Do not add nested `AGENTS.md` files.
- Do not commit, switch branches, rewrite history or perform destructive Git
  operations without explicit user instruction.
- Coordinate public API, tensor, probability, physical, schema, dependency and
  workflow changes with their governing docs and the smallest relevant check.
  Pure structural changes can use a recorded before/after comparison. Add a
  test only for an independent validation gap; do not mirror implementation or
  preserve profile-count/string assertions for coverage.
- Never silently change a physical convention. State old/new definitions and
  use a small-system reference when practical. Preserve generated research data.

## Package boundary

- MIT-licensed Julia 1.11 package, `0.6.0-DEV`, unregistered. Maintain 0.x
  semantic versioning; breaking APIs or schemas require a deliberate minor bump.
  Do not register, bundle a release or duplicate package source without a decision.
- Root dependencies: registered FiniteMPS compatible 1.x from 1.8, registered
  FiniteLattices 2.2, JLD2, Random, Dates, Printf, LinearAlgebra and SHA. Access TensorKit through
  `FiniteMPS.TensorKit`. Do not pin FiniteMPS to a path, URL, branch or commit.
- `src/StoMPS.jl` includes supported Born kernels, continuation, estimators,
  Hubbard geometry/model, mean-N thermal flow, scalar SIR/LOO/HT logic, explicit
  Hubbard sample policy, TDVP diagnostics and structured tabular logging.
- Shared orchestration and output live in `src/driver` and `src/io`. JLD2 is a
  direct dependency; MPI is optional through `ext/StoMPSMPIExt.jl`. Applications
  own parameters and MPI initialization, never duplicate runtime code or include
  another example. Scheduler configuration remains local.
- FiniteMPS owns tensor containers, contraction, SETTN and public `TDVPSweep1!`.
  Call its public interfaces from StoMPS-owned methods and keep upstream
  namespaces unchanged. Use StoMPS names for its own evolution interfaces.
  Driver step callbacks explicitly seed the default RNG; public stochastic
  kernels take an `rng::AbstractRNG`.
- Supported Born laws: rank-one local Z2 sectors; two whole charge/spin SU2SU2
  projectors; three U1SU2 occupation projectors, both exact fixed-N conditional
  and unconditioned. No generic non-Abelian, magnetic-sublevel, fusion-tree or
  multiplicity-splitting sampler without a derivation and focused validation.
- Public SU2SU2 symbols and files use `su2su2`, never legacy `so4` aliases.
  The mathematical SO(4) description may remain in physics prose.

## Tracked layout and ownership

```text
src/
test/
examples/SU2SU2_Hubbard_square_lattice_half_filling/
examples/U1SU2_Hubbard_chain_dynamics/
examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/
examples/U1SU2_Hubbard_square_lattice_grand_canonical_ensemble/
examples/U1SU2_Hubbard_triangular_lattice_canonical_ensemble/
docs/
```

- Each application has independent `config.jl`, `main.jl` and `Project.toml`.
  Source is shared only through StoMPS. Canonical apps select the shared
  ground-state DMRG; the chain selects JP dynamics. Keep their model conventions
  identical to their finite-temperature calculation.
- Only `AGENTS.md` and `README.md` are tracked root Markdown. `STRUCT.md`, plans
  and research notes are ignored. Do not replace them with tracked workflow docs.
- Z2 example, old applications, `studies/`, `ops/`, `dev/`, legacy references and
  private server support are ignored local material. Core Z2 code/tests remain
  supported. No clean-clone dependency on ignored files or personal paths.
- Generated manifests, LocalPreferences, config.local.jl, samples, checkpoints,
  figures, logs, reports and provenance remain ignored or externally archived.
  Never force generated research data into Git.

## Physical and probability contracts

- SU2SU2 factor order is charge then spin; local space is
  `su2su2_born_space()`. Whole eta/spin projectors never resolve magnetic numbers.
  Half-filled Hamiltonian includes `-mu*N`, `mu=U/2`, no constant shift; reported
  energy includes this chemical-potential term. A chain is W=1, transverse OBC.
- U1SU2 uses `u1su2_born_space()`, public order empty/single/double,
  `q=n-1`, carrier dimensions `[1,2,1]`. Canonical samples have exactly Ne;
  grand-canonical samples retain their sampled N and are never projected to Ne.
- Geometry has independent OBC/PBC axes; undirected nearest-neighbor bonds
  are unique, including length-two periodic directions. Sublattice signs come
  from physical coordinates. SU2SU2 requires a bipartite lattice. JP dynamics
  requires an open chain, centered P and J=+i[H,P].
- Triangular U1SU2 geometry is explicit; the default model remains square.
  Its bulk spin observable includes only bonds with both endpoints in bulk:
  `C_ss_tot=4*<sum_bonds S_i dot S_j>/(9*Nbulk*pbar_bulk^2)`, where
  `pbar_bulk=(<Nbulk_operator>-2*<Dbulk_operator>)/Nbulk`. Form the ratio after
  ensemble averaging. Keep `3*Nbulk`, not the number of retained bonds, as the
  spatial normalizer. Requested temperatures must match the saved beta grid.
- SETTN initializes the first beta; purification includes beta0; Born occurs
  there; subsequent cooling uses `-delta_beta/2`. All finite-T TDVP uses one-site
  CBE, default RSVD and D_cbe=D+floor(D/4). DMRG is a separate two-site reference.
- U1SU2 purification fixes mean Ne through fugacity flow. Canonical Born then
  conditions on Ne. Cross-N cooling refits lambda=beta*mu to mean Ne at every
  saved beta and once after a frozen-lambda SIR event. Missing support fails;
  never clamp mu, change Ns, switch ensembles or silently top up proposals.
- Initial Born weights are continuation norms only. Fresh defensive mixture
  branches start with log(p0)-log(Q), with Q=(1-alpha)p0+alpha*Pr for the whole
  configuration. Compression norm is diagnostic, not an importance factor.
- Physical H-only norm h accumulates independently and is never SIR-reset.
  Cross-N estimator log weight is a+(lambda-lambda0)N; TDVP adds 2log(norm) to
  a and h. SIR resets a=C-(lambda-lambda0)N, leaving h intact.
- Shared empirical proposal draws one weighted MPS member, then its whole
  configuration law, and scores the resulting configuration against all
  members. Marginalization equals the sequential posterior-RDM mixture law.
  It changes RNG consumption from historical controllers; no bitwise legacy
  trajectory compatibility is promised. MPI never averages tensors.
- Default SIR uses virtual estimator eta < 0.5, multiplicity LOO, alpha<=0.999,
  and iid pool draws until Ns distinct candidate IDs or 3Ns attempts. Repeats
  retain multiplicities; pools with fewer than Ns IDs retain all input weights
  and multiplicities without drawing. One event per beta is accepted after
  validation. SU2SU2 scheme choices remain explicit; U1SU2 supports
  with-replacement SIR/off.
  Detailed mixture, hybrid, WOR and HT contracts are in design/resampling.md.
- `Ns` is the initial count and representative target, `Nv=sum(m)` the actual
  virtual count, `K` stored representatives, `Np` distinct complete
  stored sector configurations, for every sampled beta. Keep all four labels
  distinct, including pool-ID WOR.
  Normalize estimator ESS and logZ with Nv; normalize resource eta_p with Ns.
  Legacy fixed-population modes retain Nv=Ns. The 0.6 schema rejects older snapshots.

## Persistence and reproducibility

- Two configured seeds: thermal and master. Derive operation/event/immutable
  lineage/index streams, never rank streams. Default-RNG use by FiniteMPS must
  be visibly seeded at driver boundaries and covered by rank/restart checks.
- Store config, geometry, energy convention, beta grids, sector order, cutoffs,
  seeds, status and numerical diagnostics. State total vs per-site energy.
- MPI communicates scalars and configurations only. Rank 0 prepares the beta0
  MPO on shared disk; all ranks load. Stochastic runs restart from beta0.
- Validate beta0 preparation identity and dynamics source identity on loading.
  Reusing a compatible beta0 parent may change stochastic controls and MPI size.
  No intermediate thermal/cooling tensor restore or scalar-only replay recovery.
- Save numerical history every 3600 seconds at the next safe boundary and at
  stage completion/failure. Only beta0 stores a thermal MPO. Only explicitly
  requested dynamics betas store sample MPS; DMRG stores one final MPS. Never
  put particles, environments or tensors inside a numerical result payload.
  Dynamics manifests publish after all sample files complete. Preserve old data.
- Match StoPEPS table formatting: root owns `log/run.log`; ranks own
  `log/samples/rank_XXX_sample_XXX.log`. Immutable sample indices are separate
  from mutable pool IDs and RNG lineage. No PEPS contraction columns.
- Output must be external or ignored. Runtime disk environments may use only
  explicit node-local scratch configured by launcher/TMPDIR, not persistent
  research storage by an implicit fallback. Server ABI/scheduler preferences
  are local. A system-MPI launch must validate ABI and a two-rank smoke first.

## Validation and review

- Run focused deterministic Born/normalization/weight identities, small-lattice
  ED and representative real workflows. Production benchmarks complement,
  never replace, a small independent reference.
- Retain independent numerical coverage. Shared tests cover Born/weight
  identities, ED, DMRG and JP; minimal workflows cover consecutive events,
  MPI n=1/2, beta0 reuse and numerical-only periodic persistence.
- Do not add profile-count, string-layout or exhaustive defensive-case tests.
  Keep useful ED benchmarks; archive obsolete one-off scans and server scripts.
- DMRG sweep convergence is not convergence in bond dimension or spin sector.
- Do not report completed data with nonfinite observables, wrong sectors,
  unsupported N, failed probability normalization or uncommitted history.
- Keep README useful from a clean clone; stable derivations belong in docs,
  personal deployment notes in ignored ops. Update actual API consumers when
  changing exports; expose current APIs directly and omit legacy aliases.
