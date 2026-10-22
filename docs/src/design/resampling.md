# Annealed Born resampling

Appendix B of [arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1)
constructs an improved proposal from the cooled ensemble to restore sampling
efficiency. With sample weights ``w_k``, Eq. (21) defines

```math
\eta=\frac{(\sum_{k=1}^{N_{\rm s}}w_k)^2}
{N_{\rm s}\sum_{k=1}^{N_{\rm s}}w_k^2}.
```

The empirical proposal in Eq. (23) is

```math
P_\beta(\alpha;N_{\rm s})=
\frac{\sum_k w_k\left|\langle\alpha|\psi_k(\beta)\rangle\right|^2}
{\sum_k w_k}.
```

It is combined with the initial Born distribution:

```math
P_\beta^{\rm mix}(\alpha)=
\lambda P_\beta(\alpha;N_{\rm s})+(1-\lambda)p_{\beta_0}(\alpha).
```

The mixture proportion is chosen by minimizing the second moment of the
resampled weights. An auxiliary configuration drawn from this proposal is
projected onto the purified supervector at ``\beta_0``. Its initial squared-norm
weight is ``p_{\beta_0}(\alpha)/P_\beta^{\rm mix}(\alpha)``. Imaginary-time
evolution then carries it to the current temperature, and the regenerated
samples are combined with the parent ensemble.

Applications enable resampling with `resample=true` and set its sampling
efficiency threshold with `eta_threshold`. The
[parameter reference](../guide/parameters.md) describes the calculation
controls, and [results](../guide/results.md) describes the saved output.
