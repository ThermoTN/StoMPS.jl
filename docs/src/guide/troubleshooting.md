# Troubleshooting

Start with the [2×2 serial example](quickstart.md), then use the same
configuration in MPI. Keep the original error and the stage's last committed
result. Completion flags identify finished stages.

| Error or symptom | Cause and checks | Action |
|---|---|---|
| Launcher lookup | Locate mpiexecjl in the installed depot bin directory. | Use its absolute path, usually `$HOME/.julia/bin/mpiexecjl`, or add that directory to PATH. With a custom depot, inspect `DEPOT_PATH[1]` in Julia. |
| Launcher installer says the file exists | `MPI.install_mpiexecjl()` has already been run. | Use the installed launcher; `force=true` replaces the wrapper. |
| Package resolution | Check the selected `--project` and its dependencies. | Instantiate the specific application project for its `main.jl`; use the root project for the serial library example. |
| MPI initialization/launch fails or hangs | Launcher and MPI backend may differ, or runtime resources are unavailable. | In the application environment inspect `using MPI; MPI.versioninfo()`. Try one process, then two. For custom system MPI compare launcher/library implementations and ABI; keep machine preferences local. |
| `MPI size mismatch` / `Julia thread count mismatch` | An `expected_*` check disagrees with `-n` or `-t`. | Set process/thread counts in the launcher and match the expected_* checks to them. |
| Disk environment I/O | Check TMPDIR availability, write access and free space on every node. | Set explicit node-local scratch before Julia starts. For the small example use environment_disk=false. Keep MPI output_dir shared. |
| `output exists` / unmatched `run.log` | Existing output belongs to another run. | Choose a fresh output_dir. Use overwrite_output=true to replace that output deliberately. |
| Result format check | Read the format identifier and producing package version. | Load stomps-0.6-v1 files with this package; use a matching environment for other formats. |
| `incompatible beta0 checkpoint` | Physical/thermal controls, seed or preparation grid changed. | Compare with the saved preparation fields; choose a new output directory and prepare a new parent. Stochastic Ns, master seed and MPI size are controls for each cooling run. |
| `incompatible dynamics source` / sample disagrees with manifest | The thermal/stochastic source identity or source population differs. | Use the source config and saved beta manifests with their MPS files. Independent dynamics may change time and Fourier controls. |
| `target is outside the strict sampled support` | The target mean particle number lies at the edge of, or beyond, the sampled range. | Read the reported particle-number sectors and failed beta. Choose Ns and proposal settings that sample both sides of Ne. |
| `particle-number flow failed` | The mean particle-number residual exceeds tolerance. | Inspect residuals; check thermal D, beta spacing and iteration controls against a small reference. |
| `beta0 must be a purification endpoint` / requested beta missing | beta0 or a measurement beta is inconsistent with the grid. | Include beta0 after the first beta; put dynamics/forced betas on the sampled grid. For triangular T, include beta=1/T explicitly. |
| `JP dynamics requires ... open chain` | Dynamics requested for the wrong geometry or ensemble. | Use canonical U1SU2 with W=1 and both OBC; otherwise disable dynamics_betas and ground_state_dynamics. |
| `tmax ... multiple of dt` / transform cutoff error | Time controls disagree. The cutoff error mentions ED because the transform is shared. | Set tmax to an integer multiple of dt and window_cutoff to the endpoint; transforms need at least two time points. |
| `incompatible Ne and total spin` / DMRG sweep limit | Invalid spin parity or unconverged final-D sweeps. | Check 2S parity and physical range, inspect saved sweep energies, then adjust the D/sweep schedule. Sweep success still needs D and spin-sector checks. |
| Long delay before the first numerical row | Julia compilation and model construction can dominate a small run. | Measure whole-process wall time separately. Compare a second run in a fresh output directory and inspect its evolution-step times. |
| `lnZ_error=0` or `eta_error=0` after SIR | Conditional errors can vanish for equal weights. | Assess independent-run dispersion and numerical convergence; see the uncertainty guide. |

For a two-rank MPI smoke check, after selecting the application project:

```sh
APP=examples/U1SU2_Hubbard_square_lattice_canonical_ensemble
"$HOME/.julia/bin/mpiexecjl" --project="$APP" -n 2 julia --project="$APP" -e '
using MPI
MPI.Init()
comm = MPI.COMM_WORLD
@assert MPI.Allreduce(1, +, comm) == 2
MPI.Comm_rank(comm) == 0 && MPI.versioninfo()
MPI.Finalize()
'
```

This checks the launcher and collectives; follow it
with the 2×2 application before a production run.
