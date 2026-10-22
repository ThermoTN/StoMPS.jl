# Finite-temperature evolution and weights

Imaginary-time evolution follows Eqs. (5)–(7) and App. C of
[arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1).
A normalized Born sample at ``\beta_0`` evolves as

```math
|\tilde\psi_k(\beta)\rangle=
e^{-(\beta-\beta_0)H/2}|\psi_k(\beta_0)\rangle,
\qquad w_k(\beta)=\|\tilde\psi_k(\beta)\|^2.
```

Its normalized expectation values enter the thermal average with these weights:

```math
\langle O\rangle_\beta\simeq
\frac{\sum_k w_k\langle\psi_k(\beta)|O|\psi_k(\beta)\rangle}
{\sum_k w_k}.
```

For direct continuation of the Born ensemble, the partition function obeys

```math
\ln Z(\beta)=\ln Z(\beta_0)+
\ln\mathbb E[\|\tilde\psi_k(\beta)\|^2].
```

`continue_born_sample!` stores the accumulated squared-norm logarithm in
`relative_log_weight` and normalizes the MPS after each evolution step.
`logweighted_mean` and `estimate_logpartition` evaluate the corresponding
estimators. [Annealed Born resampling](resampling.md) describes the proposal
weights for regenerated samples.

## Grand-canonical weights

For U1SU2 grand-canonical sampling, each sample retains its particle number
`N_k` during cooling under the physical Hamiltonian `H`. The chemical potential
enters the estimator through

```math
\lambda=\beta\mu,\qquad \lambda_0=\beta_0\mu_0,\qquad
w_k(\beta,\mu)=\exp\!\left[a_k(\beta)+(\lambda-\lambda_0)N_k\right].
```

For initial Born samples, `a_k` starts at zero and accumulates the logarithm
of the squared evolution norm. Resampling updates this base weight according
to the proposal and selection rules. The driver fits `lambda` after sampling
and at every saved beta, including `beta0`, so that
`sum(m_k*w_k*N_k)/sum(m_k*w_k)=Ne`; `m_k` is the stored multiplicity.
It refits after each resampling event as well. The sampled particle numbers
must bracket the target mean.

See [symmetric Born sampling](born_sampling.md#Symmetric-Born-sampling)
for the local symmetry bases and the canonical and grand-canonical treatments.

## Calculation controls

`beta_list` is a strictly increasing positive vector containing `beta0`.
SETTN initializes its first point, and tanTRG prepares the purified supervector
through `beta0`. Born sampling occurs there; subsequent MPS evolution uses
steps `-delta_beta/2` with one-site TDVP and controlled bond expansion.

SETTN, purification and samples have separate bond cutoffs. The application
config sets truncation tolerances, Krylov controls, the CBE expansion ratio
and random seeds. SU2SU2 uses the half-filled chemical potential `mu=U/2`.
U1SU2 preparation adjusts the chemical potential to the requested mean `Ne`.
Canonical Born sampling then fixes every sample to `Ne`; grand-canonical
sampling spans particle-number sectors and fits their weighted mean.

See [parameters](../guide/parameters.md) for units and settings,
[physical conventions](physical_conventions.md) for energy definitions, and
[results](../guide/results.md) for saved observables.
