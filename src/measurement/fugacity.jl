struct NumberSupportError <: Exception
    message::String
    sectors::Vector{Int}
    target::Float64
end
Base.showerror(io::IO, e::NumberSupportError) = print(io, e.message,
    "; sampled N sectors=", e.sectors, ", target=", e.target)

log_weights(base, number, lambda, lambda0) = base .+ (lambda-lambda0) .* number
physical_log_norm(norm_H, number, lambda, lambda0) =
    norm_H .+ (lambda-lambda0) .* number
sir_weight_base(common, number, lambda, lambda0) = common .- (lambda-lambda0) .* number

function number_distribution(base, number, multiplicity, lambda, lambda0)
    length(base) == length(number) == length(multiplicity) > 0 ||
        throw(DimensionMismatch("particle scalar lengths disagree"))
    all(isfinite, base) && all(>(0), multiplicity) && isfinite(lambda) && isfinite(lambda0) ||
        throw(ArgumentError("invalid weights, multiplicities or fugacity"))
    sectors = sort!(unique(Int.(number)))
    mass = [logsumexp(Float64[base[j] + log(multiplicity[j])
        for j in eachindex(base) if number[j] == n]) for n in sectors]
    function at(value)
        # Center N for exponentiation; restore the constant for absolute logXi.
        center = sectors[cld(length(sectors), 2)]
        logits = mass .+ (value-lambda0) .* (sectors .- center)
        logtotal = logsumexp(logits)
        probabilities = exp.(logits .- logtotal)
        mean_N = sum(probabilities .* sectors)
        variance_N = sum(probabilities .* (sectors .- mean_N).^2)
        return (; sectors, probabilities, mean_N, variance_N,
            log_total_weight=logtotal + (value-lambda0)*center)
    end
    return at
end

function fit_lambda(base, number, multiplicity, target, lambda0;
        initial=lambda0, tolerance=1e-8, maximum_iterations=128)
    isfinite(target) && isfinite(initial) && tolerance > 0 && maximum_iterations > 0 ||
        throw(ArgumentError("invalid fugacity solver arguments"))
    at = number_distribution(base, number, multiplicity, initial, lambda0)
    sectors = at(initial).sectors
    first(sectors) < target < last(sectors) || throw(NumberSupportError(
        "target is outside the strict sampled support; chemical potential is not identifiable",
        sectors, Float64(target)))
    left = right = Float64(initial)
    step = 1.0
    bracketed = false
    for _ in 1:64
        lo, hi = at(left).mean_N, at(right).mean_N
        if lo <= target <= hi
            bracketed = true
            break
        end
        lo > target && (left -= step)
        hi < target && (right += step)
        step *= 2
    end
    bracketed || throw(NumberSupportError("could not bracket fugacity", sectors, Float64(target)))
    value = clamp(Float64(initial), left, right)
    for iteration in 1:maximum_iterations
        distribution = at(value)
        residual = distribution.mean_N - target
        if abs(residual) <= tolerance && distribution.variance_N > eps(Float64)
            return merge(distribution, (; lambda=value, residual, iterations=iteration,
                bracket=(left,right)))
        end
        residual < 0 ? (left=value) : (right=value)
        candidate = distribution.variance_N > eps(Float64) ?
            value - residual/distribution.variance_N : NaN
        value = isfinite(candidate) && left < candidate < right ? candidate : (left+right)/2
    end
    throw(NumberSupportError("fugacity solver did not converge or variance is unresolved",
        sectors, Float64(target)))
end

function fitted_statistics(values, number, log_weight, multiplicity; beta=1.0)
    length(values) == length(number) == length(log_weight) == length(multiplicity) ||
        throw(DimensionMismatch("observable scalar lengths disagree"))
    all(isfinite, values) || throw(ArgumentError("nonfinite observable"))
    u = exp.(log_weight .- maximum(log_weight))
    u ./= sum(multiplicity .* u)
    mass = multiplicity .* u
    mean = sum(mass .* values)
    mean_N = sum(mass .* number)
    dn = number .- mean_N
    variance_N = sum(mass .* dn.^2)
    variance_N > eps(Float64) || throw(NumberSupportError(
        "variance is unresolved for fitted-mean errors", sort!(unique(number)), mean_N))
    covariance = sum(mass .* (values .- mean) .* dn)
    influence = values .- mean .- (covariance/variance_N) .* dn
    square = multiplicity .* u.^2
    correction = max(1-sum(square), eps(Float64))
    return (; mean, errorbar=sqrt(sum(square .* influence.^2)/correction),
        mu_errorbar=sqrt(sum(square .* dn.^2)/correction)/variance_N/beta,
        effective_sample_size=inv(sum(square)))
end
