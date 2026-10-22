# Library tutorials

Run these snippets in a Julia REPL started at the repository root with
`--project=.`, or save the script at the repository root.
The root environment supports serial library use. Public symbols
are listed in the [API reference](../api.md).

## A complete serial calculation

```julia
using StoMPS
include("examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/config.quickstart.jl")
c = merge(config(), (;output_dir=nothing))
result = run_application(c)
@assert result["metadata"]["complete"]
for row in result["observables"]["curve"]
    println((row.beta, row.energy_per_site, row.stderr_per_site))
end
```

`output_dir=nothing` keeps this small run in memory. The same calculation can be
split at the public preparation/sampling boundary:

```julia
parent = run_purification(c)
result = run_stochastic(c; checkpoint=parent)
```

The parent contains the normalized beta0 MPO, thermal normalization and identity.
With MPI, pass `comm=MPI.COMM_WORLD` after the application initializes MPI and
provide a shared `output_dir`; rank zero returns the numerical result.

## Draw a Born sample

Prepare a parent, choose an explicit RNG and inspect the conditional sample:

```julia
using StoMPS, Random
include("examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/config.quickstart.jl")
c = merge(config(), (;output_dir=nothing))
parent = run_purification(c)
draw = u1su2_fixed_n_purification_born_sample(parent.state, c.Ne;
    rng=MersenneTwister(42))
@assert sum(draw.configuration .- 1) == c.Ne
println((configuration=draw.configuration, log_probability=draw.log_probability))
model = hubbard_model(c)
energy_total = hubbard_expectation(draw.state, model.H; disk=false)
println(energy_total / (c.L*c.W))
```

Public configuration indices 1/2/3 mean empty/single/double, with carrier
dimensions 1/2/1. The single sector retains the full spin doublet. The direct
kernel uses its default exact projection (`notrunc()`); `hubbard_born`
applies the application's Born compression controls. The returned state's
energy is one sample's expectation value. Ensemble estimates combine sample
expectations using the squared-norm weights in Eqs. (5)–(7) of the paper.
At beta0 these weights are equal; the Born probability is represented by the
sample frequency, as in Eq. (4). See [Born sampling](../design/born_sampling.md).

For grand-canonical calculations, `u1su2_purification_born_sample(parent.state;
rng)` uses the same symmetry basis and draws across particle-number sectors.
The shared driver, `run_stochastic` with `ensemble=:grand_canonical`, then
fits the chemical potential to the target weighted mean `Ne` at each saved
beta. A direct call to the Born kernel only draws a sample; the ensemble fit
is performed by the driver.

## Load and analyze saved output

After the file-writing quick start, use the
[result-reading example](results.md) to load the saved
curve and check completion. The same guide covers fields, JP transforms and
reference files. Prefer saved ensemble fields for standard analysis;
reconstructing averages from scalars requires both log weights and multiplicities.
