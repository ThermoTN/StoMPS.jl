# Convergence and uncertainty

A completed run has passed its internal consistency checks. Establish numerical
accuracy for the observable and temperature range of interest separately.
Keep geometry, ensemble, energy convention and measurement normalization fixed
when comparing runs, and save each run in its own directory.

## Recommended sequence

1. **Check a small independent reference.** On a two-site Hubbard chain compare
   with ED or the existing deterministic Born enumeration. Canonical thermal
   results use fixed-particle-number ED including all spin sectors as their reference.
   SU2SU2 energies include the chemical-potential term. The package tests cover these conventions.
2. **Converge the thermal parent and sample representation.** Increase SETTN and
   purification bond cutoffs first, then the sample bond cutoff (for example D,
   2D, 4D until changes meet the target accuracy). Tighten truncation controls
   and check both observables and the reported maximum discarded weight.
   Changed thermal controls require a new compatible parent. DMRG additionally
   needs comparisons across final D and plausible total-spin sectors.
3. **Converge evolution grids.** Refine the preparation and cooling beta steps
   while retaining the original measurement endpoints. Check the starting beta
   and SETTN order if initialization error matters. For dynamics halve dt at
   fixed tmax, then tighten Krylov tolerance/dimension as needed. Repeat until
   shifts are below the desired physical accuracy and resolved statistical error.
4. **Measure sampling uncertainty with independent runs.** At fixed numerical
   settings use several master seeds, then increase Ns. Compare run means and
   their dispersion; a plateau as Ns grows can reveal numerical bias or poor
   population diversity. Start with 4–8 runs as a diagnostic, increasing the
   number when error estimates remain unstable. Reuse a parent to isolate
   stochastic variation; vary thermal_seed in separate parent preparations to
   assess randomized thermal compression as well.
5. **Check resampling effects.** Compare enabled/disabled resampling with matched
   numerical settings over independent seeds. Inspect sample count, eta
   and event histories. Independent runs measure the resulting uncertainty.
6. **Converge optical analysis.** Increase tmax at fixed dt with
   window_cutoff=tmax; vary broadening deliberately. A denser frequency_list
   samples the finite-time transform more finely. The time interval and window
   determine its physical resolution (paper Sec. V). Compare ED using the same
   window and broadening.

Grid refinement also changes the number of saved measurement and resampling
opportunities. First compare grids with resampling off to isolate integration
error, then check the desired adaptive workflow. Matching seeds helps comparisons
and adaptive decisions determine the resulting random trajectories.

## What the reported errors cover

`stderr_per_site` estimates sampling uncertainty using the saved multiplicities
and, for grand-canonical cooling, the fitted-mean estimator. Assess finite-D, TDVP,
SETTN, finite-size and finite-time-window effects through convergence studies.

Multiplicity records repeated selections of a sample. For R independent,
identically configured runs with estimates E_r, report the mean
and its across-run standard error

```math
\bar E=\frac1R\sum_r E_r,\qquad
\mathrm{SE}(\bar E)=\sqrt{\frac{\sum_r(E_r-\bar E)^2}{R(R-1)}}.
```

The across-run standard deviation describes the spread of individual runs.
Dividing by sqrt(R) gives the standard error of their mean. Report statistical
error and measured convergence shifts separately.

## Checks already available

The test suite contains `test/hubbard.jl` and `test/triangular_hubbard.jl`.
Run both from the repository root:

```sh
julia --project=. -t2 test/runtests.jl
```

`hubbard.jl` checks square/chain geometry, fermionic hopping conventions,
SU2SU2 and canonical U1SU2 thermal preparation and Born continuation against
independent two-site ED. Tensor evolution uses bond cutoffs of 16 and three
inverse-temperature points; the Born configurations are enumerated exactly.
It also checks multiplicity, fugacity fitting and without-replacement
resampling weight identities.

`triangular_hubbard.jl` checks triangular geometry and physical spin operators,
then compares four-site canonical calculations at Ne=3 and Ne=5 with an
independent occupation-basis ED reference. It checks energy, logZ, particle
number, double occupancy and normalized spin correlations, using purification
bond dimension 256 and sample bond dimension 64. Each test file contains its
own ED reference implementation.

Automatic CI runs both files on Julia 1.11 with registered dependencies and
two Julia threads. The test job has a 120-minute limit, including dependency
setup (30 minutes) and tests (80 minutes). Documentation builds have a separate
30-minute limit. New pushes cancel older runs on the same ref.
