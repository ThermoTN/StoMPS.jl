function _validate_multiplicity(log_weight, multiplicity, number_samples)
    length(log_weight) == length(multiplicity) || throw(DimensionMismatch(
        "log_weight and multiplicity lengths differ",
    ))
    !isempty(log_weight) && all(isfinite, log_weight) || throw(ArgumentError(
        "representative weights must be finite and nonempty",
    ))
    all(>(0), multiplicity) || throw(ArgumentError("multiplicities must be positive"))
    sum(multiplicity) == number_samples || throw(ArgumentError(
        "multiplicities must sum to number_samples",
    ))
end

function multiplicity_metrics(log_weight, multiplicity, number_samples)
    _validate_multiplicity(log_weight, multiplicity, number_samples)
    shifted = Float64.(log_weight) .- maximum(log_weight)
    logm = log.(Float64.(multiplicity))
    log_mass = shifted .+ logm
    eta = exp(2logsumexp(log_mass) -
        logsumexp(2 .* shifted .+ logm) - log(number_samples))
    eta_p = exp(2logsumexp(log_mass) -
        logsumexp(2 .* log_mass) - log(number_samples))
    eta_unique = exp(2logsumexp(log_mass) -
        logsumexp(2 .* log_mass) - log(length(log_weight)))
    return (
        eta=clamp(eta, inv(Float64(number_samples)), 1.0),
        eta_p=clamp(eta_p, inv(Float64(number_samples)),
            length(log_weight) / Float64(number_samples)),
        eta_unique=clamp(eta_unique, inv(Float64(length(log_weight))), 1.0),
        number_representatives=length(log_weight),
        unique_fraction=length(log_weight) / Float64(number_samples),
    )
end

function multiplicity_weighted_mean(values, log_weight, multiplicity, number_samples)
    _validate_multiplicity(log_weight, multiplicity, number_samples)
    length(values) == length(log_weight) || throw(DimensionMismatch(
        "values and representative weights differ in length",
    ))
    return logweighted_mean(values,
        Float64.(log_weight) .+ log.(Float64.(multiplicity)))
end

function multiplicity_weighted_statistics(
    values,
    log_weight,
    multiplicity,
    number_samples,
)
    _validate_multiplicity(log_weight, multiplicity, number_samples)
    length(values) == length(log_weight) || throw(DimensionMismatch(
        "values and representative weights differ in length",
    ))
    all(isfinite, values) || throw(ArgumentError("observable values must be finite"))
    shifted = Float64.(log_weight) .- maximum(log_weight)
    unnormalized_per_copy = exp.(shifted)
    normalization = sum(Float64.(multiplicity) .* unnormalized_per_copy)
    normalized_per_copy = unnormalized_per_copy ./ normalization
    mean = sum(Float64.(multiplicity) .* normalized_per_copy .* Float64.(values))
    sum_square_weight = sum(
        Float64.(multiplicity) .* normalized_per_copy .^ 2,
    )
    numerator = sum(
        Float64.(multiplicity) .* normalized_per_copy .^ 2 .*
        (Float64.(values) .- mean) .^ 2,
    )
    denominator = max(1 - sum_square_weight, eps(Float64))
    errorbar = sqrt(max(numerator, 0.0) / denominator)
    return (
        mean,
        errorbar,
        effective_sample_size=inv(sum_square_weight),
        sum_square_weight,
    )
end
