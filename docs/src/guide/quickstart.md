# Setup and launch

Use Julia 1.11 (`julia --version`). With Juliaup, select an installed 1.11
channel using `julia +1.11` in place of `julia`. All commands below assume the
repository root as the working directory.

## Serial quick start

The root project provides StoMPS and its numerical dependencies for serial use.
The complete configuration for 2×2 open-square Hubbard model  lives at
`examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/config.quickstart.jl`.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -t2 -e '
using StoMPS
include("examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/config.quickstart.jl")
run_application(config())
'
```

| Setting | Demonstration value |
|---|---|
| Model | L=2, W=2, OBC on both axes, t=1, U=8, canonical Ne=4 (half filling) |
| Temperature grid | beta=[1/128, 1/16, 1/8, 1/4], beta0=1/16 |
| Samples and bond dimensions | Ns=32; SETTN, purification and sample D=32 |
| Calculation | Thermal parent, Born sampling, cooling, local density and double occupancy |
| Optional stages | Resampling, JP dynamics, DMRG and ED disabled |
| Resources | One process, two Julia threads, one BLAS thread, memory environments |

The output directory (ignored in the git history tree) is
`examples/U1SU2_Hubbard_square_lattice_canonical_ensemble/output/quickstart_2x2`.
Expect `log/run.log`, sample logs, `results/checkpoint.jld2` and
`results/result.jld2`; their contents are listed in
[runtime](../design/runtime.md). Run the
[result-reading example](results.md) to check the computation, and a standard output includes: 1. a completion message;
2. three sampled temperatures; 3. the corresponding energies (without `Inf` or `Nan`); 4. exactly four particles per sample.

Installation and first-use compilation can dominate this tiny calculation.
Measure whole-process wall time with `time julia ...`; `run.log` reports
evolution-step times separately.

As a scale reference, a complete serial 2×2 validation run took 149 seconds and
peaked at 2.30 GiB resident memory on macOS arm64, Julia 1.11.6 and FiniteMPS 1.8.0
(two Julia threads, existing installed dependencies, fresh output, October 2026).
This includes package precompilation, startup, JIT compilation, the calculation,
saved-result checks and a four-site canonical ED comparison. All samples had
Ne=4, and energies at all three sampled betas agreed with ED within two sampling
standard errors. For this quick-start demonstration of the project, Ns=32 is used, 
and this is not for guarantee of computation convergence.
Budget a few minutes and several GiB; larger runs must be measured separately.
A new Julia process can repeat substantial JIT work
even after package precompilation.

To repeat a calculation, use a new `output_dir`. See
[runtime](../design/runtime.md) for output protection and beta0 reuse.

## Application environments and MPI

Each application project adds MPI and loads this checkout via `[sources]`.
Instantiate the chosen application environment before launching its entry point.
For a small MPI demonstration, make an ignored copy of the quick-start config
and change its output directory to a fresh path shared by the participating ranks:

```sh
APP=examples/U1SU2_Hubbard_square_lattice_canonical_ensemble
cp "$APP/config.quickstart.jl" "$APP/config.local.jl"
julia --project="$APP" -e 'using Pkg; Pkg.instantiate()'
```

`main.jl` initializes MPI. A one-process MPI
run is:

```sh
julia --project="$APP" -t2 "$APP/main.jl" --config "$APP/config.local.jl"
```

For multiple processes, install the environment-aware launcher once:

```sh
julia --project="$APP" -e 'using MPI; MPI.install_mpiexecjl()'
```

Its default location is `joinpath(DEPOT_PATH[1], "bin", "mpiexecjl")`, usually
`$HOME/.julia/bin/mpiexecjl`. Use that absolute path or add its directory to
`PATH`. If already installed, skip installation. With the default depot:

```sh
export PATH="$HOME/.julia/bin:$PATH"
mpiexecjl --project="$APP" -n 2 julia --project="$APP" -t2 "$APP/main.jl" --config "$APP/config.local.jl"
```

Use another fresh output directory when comparing one and two processes.
Keep the same seeds and numerical parameters; rank zero owns the main
result. `-n` sets MPI size, `-t` Julia threads per process, and `blas_threads=1`
sets one BLAS thread per rank.

## Disk environments and production settings

The quick start uses `environment_disk=false`. For `environment_disk=true`,
set `TMPDIR` **before starting Julia** to an existing, writable node-local
scratch directory. On a local workstation, for example:

```sh
export TMPDIR="$(mktemp -d /tmp/stomps-env.XXXXXX)"
```

On a cluster, use the scheduler's node-local scratch allocation on every node.
Keep `output_dir`
on shared research storage for MPI checkpoints and results. Scratch is temporary
contraction storage and is separate from saved results.

Use the MPI backend configured for the application environment. A custom
system MPI uses a matching launcher and Julia MPI backend. Check the ABI and
run a two-rank smoke calculation before production. See
[troubleshooting](troubleshooting.md).

For a larger calculation, copy that application's own `config.jl` to
`config.local.jl`, edit the [parameters](parameters.md), and pass it
with `--config`. And please alway make sure the results are converged when varying the parameter (see [convergence workflow](validation.md) ).
Available task types and their outputs are in the
[task table](../design/runtime.md).
