"""
    logweighted_mean(values, log_weights)

Compute a weighted mean when the non-negative weights are represented by their
logarithms.  A maximum-log shift is applied before exponentiation, so a common
additive offset in `log_weights` leaves the result invariant.

For samples drawn directly from the exact beta0 Born distribution, that Born
probability is represented by sample frequency. Their weights are the squared
norms accumulated during continuation. A sample regenerated from a mixture
proposal `P` starts with log weight `log(p_beta0) - log(P)` before continuation
increments are added. See Eqs. (4)–(7) and App. B of
[arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1).
"""
function logweighted_mean(
    values::AbstractVector{<:Number},
    log_weights::AbstractVector{<:Real},
)
    _validate_logweighted_inputs(values, log_weights)
    shift = maximum(log_weights)
    weights = exp.(log_weights .- shift)
    return sum(weights .* values) / sum(weights)
end

"""
    estimate_logpartition(logZ0, relative_log_weights)

Estimate the absolute log partition function from the thermal-MPO value
`logZ0` at the Born sampling inverse temperature beta0 and the samples' continuation log
weights.  The result is

```text
logZ0 + logmeanexp(relative_log_weights).
```
"""
function estimate_logpartition(
    logZ0::Real,
    relative_log_weights::AbstractVector{<:Real},
)
    isfinite(logZ0) || throw(ArgumentError("logZ0 must be finite"))
    return logZ0 + _logmeanexp(relative_log_weights)
end

"""
    sampling_efficiency(log_weights)

Return the sampling efficiency eta of Eq. (21) in
[arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1).
For a vector of finite log weights of `N_s` samples, `N_s >= 1`, this is

```text
exp(2 * logsumexp(ell) - logsumexp(2ell) - log(N_s)).
```

The result is stable under a common additive shift of all weights and is
clamped to `[1 / N_s, 1]` within floating-point roundoff. Initial Born samples
use continuation squared-norm weights. Regenerated samples include the
mixture proposal correction `log(p_beta0) - log(P)` and subsequent squared-norm
increments.
"""
function sampling_efficiency(log_weights::AbstractVector{<:Real})
    isempty(log_weights) && throw(ArgumentError(
        "at least one log weight is required",
    ))
    all(isfinite, log_weights) || throw(ArgumentError(
        "log weights must be finite",
    ))
    count = length(log_weights)
    log_ess = 2 * logsumexp(log_weights) -
        logsumexp(2 .* log_weights) - log(count)
    efficiency = exp(log_ess)
    lower = inv(Float64(count))
    roundoff = 1e-14
    lower - roundoff <= efficiency <= 1 + roundoff || error(
        "sampling efficiency $efficiency is outside [$lower, 1]",
    )
    return clamp(efficiency, lower, 1.0)
end

function _logmeanexp(values::AbstractVector{<:Real})
    isempty(values) && throw(ArgumentError("at least one log weight is required"))
    all(isfinite, values) || throw(ArgumentError("log weights must be finite"))
    return logsumexp(values) - log(length(values))
end

function _validate_logweighted_inputs(values, log_weights)
    isempty(values) && throw(ArgumentError("at least one value is required"))
    length(values) == length(log_weights) || throw(DimensionMismatch(
        "values and log_weights must have the same length",
    ))
    all(isfinite, values) || throw(ArgumentError("values must be finite"))
    all(isfinite, log_weights) || throw(ArgumentError("log weights must be finite"))
    return nothing
end
