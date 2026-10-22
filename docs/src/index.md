# StoMPS.jl

```@meta
CurrentModule = StoMPS
```

User and method documentation for StoMPS Hubbard calculations.

1. Run the [2×2 quick start](guide/quickstart.md), then
   [read its results](guide/results.md).
2. Choose an [application](guide/applications.md) and consult the
   [parameter reference](guide/parameters.md) before increasing the system size.
3. Follow [convergence and uncertainty](guide/validation.md) to assess accuracy,
   and use [troubleshooting](guide/troubleshooting.md) for launch or calculation errors.

## Method and implementation

- [Paper and notation](design/paper.md)
- [Physical and lattice conventions](design/physical_conventions.md)
- [Tensor and symmetry conventions](design/tensor_conventions.md)
- [Born sampling](design/born_sampling.md)
- [Thermal evolution and weights](design/finite_temperature.md)
- [Execution, output and method interfaces](design/runtime.md)
- [Resampling methods](design/resampling.md)

## Library API

The [library tutorials](guide/library.md) show serial orchestration and an
independent Born draw. Use the [API reference](api.md) for details.
