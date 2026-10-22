"""
    continue_born_sample!(sample, beta_grid; rng, step!, observe,
                          normalization_tolerance=1e-10)

Continue one normalized Born sample from the first inverse temperature in
`beta_grid` to every later point.  StoMPS owns the finite-temperature weight
convention while the caller supplies the concrete evolution and observation
operations:

- `step!(sample, -delta_beta / 2, beta_index, rng)` mutates `sample` by one
  imaginary-time step;
- `observe(sample, beta, beta_index)` returns one finite scalar observable.

After every raw step, this function adds `2 * log(norm(sample))` to the relative
log weight and normalizes the sample. A step using `evolution_step!` already
normalizes and returns `log_norm_increment`; that increment is consumed once.
The initial relative log weight is exactly zero because Born probability is
already represented by sample frequency.

The weight convention follows Eqs. (4)–(7) of
[arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1).

The explicit RNG is passed to `step!`; that callback owns any default-RNG
seeding needed by FiniteMPS. The returned named tuple has fields `beta`, `observable`, and
`relative_log_weight`.  The input `sample` remains normalized at the final beta.
"""
function continue_born_sample!(
    sample::DenseMPS,
    beta_grid::AbstractVector{<:Real};
    rng::AbstractRNG,
    step!::Step,
    observe::Observe,
    normalization_tolerance::Real=1e-10,
) where {Step,Observe}
    beta = Float64.(beta_grid)
    _validate_continuation_beta_grid(beta)
    normalization_tolerance > 0 || throw(ArgumentError(
        "normalization_tolerance must be positive",
    ))

    initial_norm = norm(sample)
    isfinite(initial_norm) || throw(ArgumentError("sample norm must be finite"))
    isapprox(initial_norm, 1; atol=normalization_tolerance, rtol=0) || throw(
        ArgumentError("Born continuation requires a normalized sample"),
    )

    first_observable = observe(sample, beta[1], 1)
    _validate_observable(first_observable, 1)
    observable = Vector{typeof(first_observable)}(undef, length(beta))
    observable[1] = first_observable
    relative_log_weight = zeros(Float64, length(beta))

    for beta_index in 2:length(beta)
        delta_beta = beta[beta_index] - beta[beta_index - 1]
        diagnostics=step!(sample, -delta_beta / 2, beta_index, rng)
        increment=hasproperty(diagnostics,:log_norm_increment) ? diagnostics.log_norm_increment :
            normalize_logweight!(sample)
        relative_log_weight[beta_index] =
            relative_log_weight[beta_index - 1] + increment
        isapprox(norm(sample), 1; atol=normalization_tolerance, rtol=0) || error(
            "Born continuation failed to normalize the sample at beta " *
            "$(beta[beta_index])",
        )

        value = observe(sample, beta[beta_index], beta_index)
        _validate_observable(value, beta_index)
        observable[beta_index] = value
    end

    return (; beta, observable, relative_log_weight)
end

function _validate_continuation_beta_grid(beta)
    isempty(beta) && throw(ArgumentError("beta_grid must be nonempty"))
    all(isfinite, beta) || throw(ArgumentError("beta_grid must be finite"))
    all(value -> value >= 0, beta) || throw(ArgumentError(
        "beta_grid must contain non-negative inverse temperatures",
    ))
    all(diff(beta) .> 0) || throw(ArgumentError(
        "beta_grid must be strictly increasing",
    ))
    return nothing
end

function _validate_observable(value, beta_index)
    value isa Number || throw(ArgumentError(
        "observe must return a scalar Number at beta index $beta_index",
    ))
    isfinite(value) || throw(ArgumentError(
        "observe returned a non-finite value at beta index $beta_index",
    ))
    return nothing
end
