# Finite-temperature evolution and weights

`beta_list` is a strictly increasing positive vector. The default is

```julia
vcat(2.0 .^ vcat(collect(-15:1:-2),collect(-1.5:0.5:1)),collect(3.0:1:16))
```

SETTN initializes its first point. Purification includes beta0 (default 1),
which must appear in the vector. Born sampling occurs at beta0, and later
points use normalized MPS continuation with step -delta_beta/2. Every
finite-temperature TDVP is one-site CBE; default RSVD is on and
D_cbe=D+floor(D*cbe_ratio), cbe_ratio=0.25. SETTN, purification and samples
have separate bond cutoffs. Every tolerance, Krylov dimension, SETTN noise
schedule and seed is explicit in config.jl.

U1SU2 starts from exp(-alpha*N/2), where alpha=log((2-n)/n), n=Ne/(L*W).
This analytic choice sets the infinite-temperature density. Coupled H/N tanTRG
adjusts alpha and corrects residual N when needed. Its saved lambda is -alpha
and mu=lambda/beta.
SU2SU2 directly uses mu=U/2. Starting from normalized identity, its logZ
includes the initial Hilbert-space dimension plus twice every norm logarithm.

An initial normalized Born sample has log estimator weight a=0 and physical
norm h=0. TDVP contributes 2log(norm) to both, then the state is normalized.
No compression norm is added to either physical or importance weights.
Fresh proposal states instead start at a=log(p0)-log(Q), h=0 and replay the
same beta0-to-current cooling grid. This is the explicit-mixture importance
correction, applied exactly once.

Canonical partition normalization is
logZ_N(beta0)=logXi(beta0,lambda0)+log Pr(N=Ne)-lambda0*Ne.
SU2SU2 starts from its grand-canonical logZ. At later beta, both add
logsumexp(log(m)+a)-log(Nv), where Nv=sum(m) is the actual virtual count.
Canonical energy is H only; SU2SU2 energy is K.

For cross-N cooling, per-copy estimator log weights are
ell_j=a_j+(lambda-lambda0)N_j. Fit lambda so their multiplicity-weighted
mean equals Ne. Physical LOO norm is h_j+(lambda-lambda0)N_j, independent
of SIR resets. SIR freezes lambda, resets a_j=C-(lambda-lambda0)N_j in sampled
equal-weight branches (undersized adaptive pools retain a_j), then
refits lambda once. Eta may be below one afterward. Every saved beta meets
the requested mean-N tolerance; unsupported targets fail without committing.

The two seeds are configured by the application. FiniteMPS randomized operations
receive explicit driver default-RNG seeds derived from operation,
generation, immutable lineage and step, independently of MPI rank.
Born/resampling draws accept explicit rng objects. The empirical mixture's
member-first draw and probability law are described in
[resampling.md](resampling.md).

Results record total sample energy and per-site ensemble energy, Ns/Nv/K/Np,
mu, eta, eta_p, logZ, weight/proposal/norm provenance, local observables and
diagnostics. Scalar history includes every completed beta and is saved hourly.
See [execution and output](runtime.md) for file contents, restart behavior and
the versioned storage contract.

## Scheduled bulk spin measurements

The triangular canonical application measures energy at every sample beta.
Each requested correlation temperature T (k_B=1, in the same energy units as t)
must match exactly one saved sample beta via beta=1/T, allowing only floating
roundoff. The application does not insert beta points or interpolate data.
An empty temperature list disables correlations. Measurements occur after
any SIR event at that beta, so observables and saved weights describe the same
population. They neither draw random numbers nor alter physical norm weights.

For each stored representative k, retain q_k=<Q_B>_k, n_k=<sum_B n_i>_k and
d_k=<sum_B n_up*n_down>_k, aligned by particle ID with its scalar history.
Let p_k=(n_k-2d_k)/N_B and
`w_k=m_k*exp(a_k)/sum_j(m_j*exp(a_j))`. First form the weighted means qbar and
pbar, then report `f=A*qbar/pbar^2`, with `A=4/(9*N_B)`. Averaging the individual
ratios would change the observable. Individual representatives may have zero
singles probability; only the ensemble denominator must be positive.

The quoted first-order error propagates the joint numerator and denominator
through the centered influence

```math
g_k=A\left[\frac{q_k-\bar q}{\bar p^2}
 -\frac{2\bar q(p_k-\bar p)}{\bar p^3}\right].
```

Use the existing multiplicity-weighted particle error estimator on g_k.
This retains numerator/denominator covariance and the virtual-copy meaning of
m_k. Store the raw moments and their scalar weights for independent-run
uncertainty analysis. Per-bond scalars are
not stored. Geometry, normalization definition and temperature schedule are
part of the numerical provenance and dynamics-source identity.
