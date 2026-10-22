# Reading results

Run Julia from the repository root with `--project=.`. Numerical files use
`stomps-0.6-v1`; `load_result` checks that format. Dictionaries use string keys;
configuration and curve rows are named tuples accessed with dots.

## Energy and completion

After the [quick start](quickstart.md), run this in a Julia REPL started at the
repository root, or in a script saved there:

```julia
using StoMPS
include("examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/config.quickstart.jl")
c = config()
r = load_result(joinpath(c.output_dir, "results", "result.jld2"))
meta = r["metadata"]
@assert meta["complete"] && meta["stage_complete"]
rows = r["observables"]["curve"]
@assert [x.beta for x in rows] == sample_betas(r["parameters"])
@assert all(x -> isfinite(x.energy_per_site) && isfinite(x.stderr_per_site), rows)
@assert all(p.N == c.Ne for x in rows for p in x.scalars)
beta = [x.beta for x in rows]
temperature = 1 ./ beta
energy = [x.energy_per_site for x in rows]
stderr = [x.stderr_per_site for x in rows]
foreach(println, zip(temperature, energy, stderr))
local_density = last(rows).local_observables.number
double_occupancy = last(rows).local_observables.double
```

The assertion above is for completed **cooling**. For a purification run,
`stage == "purification"` and `stage_complete == true` establish stage success;
`complete` remains false until cooling finishes. A full purification
reference has its own file and completed flag. Other complete trajectories set
both flags; a stopped cooling run remains incomplete.
For stopped/failed runs, inspect `metadata["progress"]`,
`get(metadata, "error", nothing)` and the last committed curve row.
The committed curve records completed ensembles.

With `task=:all`, check cooling completion in `result.jld2` and requested
dynamics completion in `dynamics_result.jld2`.
The return value of `run_application` for `:all` is the cooling result.
File names and restart rules are in [runtime](../design/runtime.md).

## Common field dictionary

| Field | Type / units | Meaning |
|---|---|---|
| `parameters` | NamedTuple | Effective configuration for this stage. |
| `metadata` | Dictionary | Stage, `complete`/`stage_complete` flags, progress, run identity and provenance. |
| `metadata["provenance"].energy_convention` | String | Hamiltonian and chemical-potential convention; see [physical conventions](../design/physical_conventions.md). |
| `observables["curve"]` | Vector of NamedTuples | Committed sampling/cooling rows, including beta0. |
| `observables["purification"]` | Vector of NamedTuples | Grand-canonical thermal preparation rows. |
| `row.beta` | Float64, inverse energy | Saved inverse temperature. |
| `row.energy_per_site`, `row.stderr_per_site` | Float64, energy/site | Ensemble energy and sampling standard error. |
| `row.mean_N`, `row.mu` | Real particle count; energy or nothing | Canonical row.mu is nothing; fitted mu belongs to grand-canonical cooling or thermal preparation. |
| `row.logZ`, `row.lnZ_error` | Float64, dimensionless | Log partition function and conditional error. |
| `row.eta`, `row.eta_error` | Float64, dimensionless | Sampling efficiency and its conditional error. |
| `row.scalars` | Vector of NamedTuples | Per-sample identifiers, multiplicity, configuration, particle number, total energy and weight data. |
| `row.log_weight` | Vector of Float64 | Per-copy estimator log weights aligned with scalars; multiply mass by m when forming ensemble averages. |
| `row.local_observables` | NamedTuple or nothing | Site vectors `number`, `double`, `number_stderr` and `double_stderr`, in MPS site order. |
| `row.maximum_truncation`, `row.maximum_entropy`, `row.maximum_bond_dimension` | Real, real, integer | Evolution diagnostics; truncation is squared second-half-sweep TrunErr. |
| `row.elapsed_seconds` | Float64, seconds | Time spent on the evolution step. |
| `diagnostics["events"]` | Vector | Resampling audits. |

`nothing` usually means an optional measurement is disabled, scheduled at another beta, or a
quantity is inapplicable. `NaN` is used for undefined diagnostics, such as initial
TDVP truncation/entropy or an unavailable efficiency error. Preserve these values
when reading diagnostics. Nonfinite physical observables indicate a failed
calculation. Logs print some nothing fields as NaN; the JLD2 file retains
their types. See [uncertainty and convergence](validation.md).

## JP dynamics and references

After a chain run with saved sources and completed dynamics:

```julia
using StoMPS
d = load_result(joinpath("examples/U1SU2_Hubbard_chain_dynamics/output",
                         "results", "dynamics_result.jld2"))
@assert d["metadata"]["complete"]
c = d["parameters"]
for beta in c.dynamics_betas
    i = findfirst(==(beta), sample_betas(c))
    response = d["observables"]["beta_$(i)"]
    time, JP = response.time, response.Cgreater_JP
    println((beta, t0=JP[1], t0_check=response.t0_residual_abs))
    if response.transform !== nothing
        omega = c.frequency_list
        sigma = response.transform.conductivity
        foreach(println, zip(omega, real.(sigma), imag.(sigma)))
    end
end
```

The key uses the **sample-grid index**. JP is a
complex vector of total extensive correlations; `stderr_real` and `stderr_imag`
are time-domain sampling errors. `transform` contains kernel, damping and
complex conductivity per chain length, or is nothing when `frequency_list` is
empty. The stored per-sample trajectories support frequency-domain uncertainty
analysis that accounts for the correlations between time points sharing the
same samples. The [optical convention](../design/physical_conventions.md)
defines signs, units and size normalization.

`ground_state_result.jld2` has `observables["sweeps"]`, `observables["ground_state"]`
and optional `observables["dynamics"]`. The latter stores connected JP time data,
for subsequent conductivity analysis. The summary distinguishes
sweep convergence from bond-dimension and spin-sector convergence.
`ed_result.jld2` stores `observables["reference"]`; its `analyses` vector contains
beta, total energy, correlation, kernel and conductivity. Divide its `energy_total`
by LW before comparing with `energy_per_site`. Tensor files and the beta0 restart
procedure are listed in [runtime](../design/runtime.md).
