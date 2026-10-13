# Resampling contracts

The shared numerical implementation is `src/measurement/`,
`src/measurement/fugacity.jl` and `src/resample/`. The shared driver owns MPI
coordination and output. SU2SU2 provides the schemes described below; the
U1SU2 applications expose only with-replacement multiplicity LOO-SIR or off.

## Counts and weights

Ns is the initial Born count and the target number of stored MPS representatives.
A stored representative j has m_j positive integer virtual copies, sum m_j=Nv;
K is the number of stored representatives. Fixed-count modes have Nv=Ns. With
per-copy weight w_j and aggregate W_j=m_j*w_j:

```math
\eta={ (\sum_j W_j)^2\over N_v\sum_j W_j^2/m_j},\qquad
\eta_p={ (\sum_j W_j)^2\over N_s\sum_j W_j^2}.
```

Eta is the expanded-population estimator ESS divided by actual Nv; eta_p is
the candidate-ID ESS divided by target Ns, used for augmentation. The low-level
`multiplicity_metrics` normalizes both ESS values by its actual-count argument;
application sample statistics rescale its eta_p by Nv/Ns. Np counts different
**complete sector configurations**, including in pool-ID WOR. Np need not equal K or Ns.
Before Born sampling, Np is absent in logs. Every subsequent beta logs Ns/Nv/K/Np.

The default trigger is strict eta<0.5, resampling enabled, one event at a saved
beta and unconditional acceptance after validation. A low post-event eta
never causes another event at the same beta. An event budget limits distinct
beta events. Turning `resample=false` retains the same thermal parent, Born
seed, sample count and cooling parameters for a no-resampling control.

## Empirical proposal and defensive mixture

Pr(X) is the weighted mixture of current normalized MPS configuration laws.
Multiplicity contributes to its weights. Draw a member according to W_j,
then draw X sequentially from that member. Scoring X against every member
gives the actual Pr(X). This marginal law is identical to the
posterior-weighted RDM mixture at every site:

```math
P_r(X)=\sum_j {W_j\over\sum_k W_k}P_j(X).
```

Owners draw and score their local MPS; MPI collectives exchange scalar scores
and complete sector configurations. Member-first and posterior-RDM sampling
consume random numbers in different orders, so trajectory reproducibility
depends on the sampling algorithm as well as the seed.

Each fresh full configuration uses Q_alpha=(1-alpha)p0+alpha*Pr with one
component draw for the **whole configuration**. Alpha=1 is the empirical
endpoint; alpha=0 is baseline Born. Alpha=1 removes the baseline support
guarantee: the importance identity then applies only on the support covered
by the empirical ensemble. For a projected state, compression follows
normalization; its retained compression norm is diagnostic. Initialize a with
log(p0)-log(Q_alpha), then add all beta0-to-current norm increments. Baseline
Born draws reuse their projected states. Q_alpha uses the current empirical
ensemble, and its importance correction is applied once.

## Alpha selection

- `:fixed` uses `resampling_alpha`; alpha=1 retains the pure empirical scheme.
- `:loo_first` scores the old configurations after excluding their MPS member
  and is restricted to a single event. It minimizes empirical second-moment
  risk over [0,alpha_max]. It draws no validation samples or configuration
  histogram. Equal configurations can still be different ensemble members.
- `:loo_generalized` uses squared amplitude
  `2a+log(Q_previous)-log(p0)` for later generations, while the empirical
  proposal still uses the unmodified estimator weights. It supports repeated
  all-regenerate events and the explicit experimental eta-improvement policy.
- `:loo_multiplicity` excludes the **whole representative lineage**, uses
  empirical coefficients log(m_j)+ell_j, and risk log amplitude
  `(log(m_j)+ell_j+h_j+(lambda-lambda0)*N_j)/2` for cross-N. The last fugacity
  term is absent for SU2SU2/fixed-N. Independent h is never SIR-reset.

The common risk is the mean squared-amplitude factor divided by
`(1-alpha)+alpha*Pr_LOO/p0`, evaluated in log space. The bounded convex search
also compares both endpoints. Default alpha_max=.999, tolerance=1e-8;
single-lineage LOO falls back to alpha=0. The audit saves scores, selected
alpha, risks, pre/fresh/post weights and configurations.

## All regeneration and hybrid

`:regenerate` replaces all Ns branches with fresh mixture draws, each with its
actual proposal correction and full replayed physical norm. Repeated draws
remain separate lineages in this mode.

`:hybrid` keeps floor(eta*Ns) draws, clipped to 1:Ns-1, and regenerates the rest.
Keep ancestors are independent inverse-CDF draws of aggregate estimator mass;
each receives the mean old estimator weight. Fresh draws retain their mixture
correction and physical norm. The combined Ns population immediately estimates
observables and logZ from its accumulated estimator weights. Eta is computed
from those weights.

For non-SIR SU2SU2 calculations, `eta_weight=:physical_norm` selects
`a-log(p0)+log(Q)` for triggering; observables/logZ still use a. The default is
`:estimator`. SIR requires estimator eta.

## Default adaptive multiplicity SIR with replacement

The default `sir_draw_mode=:with_replacement_until_unique` generates
M=clamp(Ns-floor(ESS_ID),0,Ns-1), where ESS_ID=eta_p*Ns. This resource target
does not grow with the actual virtual count. Fresh draws are grouped by complete
configuration before evolution; old IDs and fresh IDs remain distinct even on
equal configurations. The pool has old masses m*w and fresh masses r*w,
and actual virtual count Nv_pool=Nv_pre+M.

If the pool has fewer than Ns candidate IDs, keep every ID, weight and input
multiplicity without consuming any selection RNG draws. Otherwise draw iid
IDs with fixed probability mass/sum(mass), stopping at Ns distinct IDs or 3Ns
attempts. A pool of exactly Ns IDs still enters this loop. Repeated hits count
toward the attempt budget and are stored as multiplicities. Reaching the target
on the last allowed draw is recorded as target_reached; a cap hit with fewer
IDs is accepted without top-up. The selected per-copy weight is
sum(pool mass)/Nv_pool, and Nv_post equals the number of attempts. An undersized
pool instead has Nv_post=Nv_pool and is not limited to 3Ns virtual copies.

Record the draw sequence, draw_count, stop_reason (target_reached, draw_limit,
pool_under_target), K_pool, Nv_pre, Nv_pool and Nv_post in the event audit.
Each selected MPS evolves once, applying its increment to all virtual copies.
All observable diagnostics and logZ use Nv=sum(m); logZ adds
logsumexp(log(m)+ell)-log(Nv) to logZ0. This preserves pool mean-weight scale.

## Fixed-count multiplicity SIR with replacement

`:sir`, `sir_draw_mode=:with_replacement` creates
M=Ns-floor(eta_p*Ns) fresh mixture draws. When M=0, only the old population
enters the pool. Equal fresh complete configurations
are grouped before evolution, using the first draw's deterministic lineage;
old and fresh candidate IDs remain distinct even if their configurations agree.
Old candidates have mass m*w; fresh candidates have mass r*w. Draw Ns virtual
slots with replacement from this combined pool. Selected slots share per-copy
weight sum(pool mass)/(Ns+M). Repeated draws of a candidate ID become one
representative with multiplicity. Subsequent TDVP evolves it once and applies
its norm increment to all virtual copies. Observables/logZ use expanded Ns.

## Configuration WOR

`:without_replacement_configuration` creates exactly M=Ns-floor(eta*Ns) fresh
draws and groups **all** old/fresh complete configurations. No adaptive top-up,
minimum pool ratio, rejection-to-uniqueness or change of Ns is allowed. The
initial iid Born ensemble is unchanged. A pool with Np<Ns fails **before fresh
regeneration**, preserving the last numerical result; equality warns and keeps
all pool configurations.

A genuinely new configuration evolves on its first fresh draw's lineage.
For an old hit, reuse the first old global ID's MPS and independent h. Its
fresh per-draw weight is exp(log(p0)-log(Q)+h), including fugacity for cross-N
where applicable; add r times this weight to the old aggregate mass. Never
multiply an old HT-corrected weight by fresh repetition count. State reuse is
an explicit numerical convention because independent truncation may differ.

For aggregate configuration masses A_i, choose capped inclusion probabilities
pi_i=min(1,A_i/tau), sum pi_i=Ns. An independently seeded random permutation
and systematic offset select Ns configurations. Each retained configuration
has multiplicity one and weight `(Ns/(Ns+M))*A_i/pi_i`. This preserves pool
mean-weight normalization and conditional weighted totals. It does not replace
physical norms or proposal probabilities with inclusion probabilities. Eta
equals eta_p and need not be one. Valid events commit once regardless of eta.

Persist raw/merged masses, configurations, grouping and first-draw maps,
cached old IDs, proposal scores and the full inclusion design. Selected
complete configurations must be distinct.

## Pool-ID WOR

`:without_replacement_pool_id` uses eta for triggering and eta_p for the fresh
count. Group fresh repeats only; old IDs remain distinct and a fresh old hit
is evolved independently. Apply the same capped/systematic/HT design to old
m*w and fresh r*w masses without merging configurations. Ns different IDs
survive with multiplicity one; their configurations can repeat. Its pool size
is an ID count K_pool, while the unified Np still means distinct configurations.

## Trigger choices

SU2SU2 retains `:threshold_each`, `:threshold_first`, `:fixed_beta`,
`:forced_and_threshold` and `:threshold_until_target`. The latter stops future
events after an event reaches eta_target; it never retries the same beta.
`forced_betas` must be explicit members of the saved beta grid. The
eta-improvement acceptance choice is confined to generalized-LOO all
regeneration; the production default always accepts valid events.

## Cross-N weights

Cross-N weights are ell=a+(lambda-lambda0)N. Physical H-only h is accumulated
independently. Freeze lambda during each SIR event; in the equal-weight branches
reset a=C-(lambda-lambda0)N (undersized adaptive pools retain their input a),
then fit lambda once to mean Ne. The new eta need not be one. A missing sampled
N support fails instead of clamping mu or projecting to fixed N.

The shared driver saves numerical history and event audits according to the
[output contract](runtime.md).
