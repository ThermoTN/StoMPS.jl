function logsumexp(values::AbstractVector{<:Real})
    isempty(values) && throw(ArgumentError("logsumexp input must be nonempty"))
    shift = maximum(values)
    shift == -Inf && return -Inf
    return shift + log(sum(exp(value - shift) for value in values))
end

function linear_mixture_log_probability(log_p0::Real, log_empirical::Real, alpha::Real)
    0 <= alpha <= 1 || throw(ArgumentError("alpha must lie in [0,1]"))
    alpha == 0 && return Float64(log_p0)
    alpha == 1 && return Float64(log_empirical)
    return logsumexp(Float64[
        log1p(-alpha) + log_p0,
        log(alpha) + log_empirical,
    ])
end


function sir_candidate_count(eta_p::Real, number_samples::Integer)
    0 < eta_p <= 1 || throw(ArgumentError("eta_p must lie in (0,1]"))
    number_samples > 1 || throw(ArgumentError("SIR requires at least two samples"))
    return clamp(number_samples - floor(Int, eta_p * number_samples),
        0, number_samples - 1)
end

function deduplicate_configurations(configurations::AbstractMatrix{<:Integer})
    number_draws, number_sites = size(configurations)
    number_draws > 0 && number_sites > 0 || throw(ArgumentError(
        "configuration matrix must be nonempty",
    ))
    by_key = Dict{Tuple,Int}()
    rows = Vector{Vector{Int}}()
    first_draw_id = Int[]
    multiplicity = Int[]
    draw_to_representative = zeros(Int, number_draws)
    for draw_id in 1:number_draws
        row = Int.(vec(configurations[draw_id, :]))
        representative = get(by_key, Tuple(row), 0)
        if representative == 0
            representative = length(rows) + 1
            by_key[Tuple(row)] = representative
            push!(rows, row)
            push!(first_draw_id, draw_id)
            push!(multiplicity, 0)
        end
        multiplicity[representative] += 1
        draw_to_representative[draw_id] = representative
    end
    unique_rows = zeros(Int, length(rows), number_sites)
    for representative in eachindex(rows)
        unique_rows[representative, :] .= rows[representative]
    end
    return (
        configurations=unique_rows,
        draw_to_representative,
        first_draw_id,
        multiplicity,
        number_draws,
        number_representatives=length(rows),
    )
end

function _risk(log_amplitude, log_ratio, alpha)
    terms = Float64[
        2log_amplitude[index] -
        linear_mixture_log_probability(0.0, log_ratio[index], alpha)
        for index in eachindex(log_amplitude)
    ]
    return logsumexp(terms) - log(length(terms))
end

function loo_alpha_selection(log_amplitude, log_p0, loo_log_empirical;
        alpha_max=0.999, tolerance=1e-8)
    length(log_amplitude) == length(log_p0) == length(loo_log_empirical) ||
        throw(DimensionMismatch("LOO arrays must have equal length"))
    isempty(log_amplitude) && throw(ArgumentError("LOO arrays must be nonempty"))
    length(log_amplitude) == 1 && return (
        selected_alpha=0.0,
        risk_at_zero=_risk([0.0], [-Inf], 0.0),
        risk_at_selected=_risk([0.0], [-Inf], 0.0),
        log_density_ratio=[-Inf],
        fallback="single-lineage",
    )
    log_ratio = Float64.(loo_log_empirical) .- Float64.(log_p0)
    amplitude = Float64.(log_amplitude) .- maximum(log_amplitude)
    left, right = 0.0, Float64(alpha_max)
    golden = (sqrt(5.0) - 1) / 2
    x1 = right - golden * (right - left)
    x2 = left + golden * (right - left)
    f1, f2 = _risk(amplitude, log_ratio, x1), _risk(amplitude, log_ratio, x2)
    while right - left > tolerance
        if f1 <= f2
            right, x2, f2 = x2, x1, f1
            x1 = right - golden * (right - left)
            f1 = _risk(amplitude, log_ratio, x1)
        else
            left, x1, f1 = x1, x2, f2
            x2 = left + golden * (right - left)
            f2 = _risk(amplitude, log_ratio, x2)
        end
    end
    candidates = Float64[0.0, alpha_max, (left + right) / 2]
    risks = [_risk(amplitude, log_ratio, alpha) for alpha in candidates]
    selected_index = argmin(risks)
    return (
        selected_alpha=candidates[selected_index],
        risk_at_zero=risks[1],
        risk_at_selected=risks[selected_index],
        log_density_ratio=log_ratio,
        fallback="none",
    )
end

function _splitmix64(value::UInt64)
    z = value + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

function deterministic_seed(master_seed::Integer, domain::UInt64,
        event_id::Integer, item_id::Integer)
    value = UInt64(master_seed) ⊻ domain ⊻ UInt64(event_id) * 0x9e3779b97f4a7c15 ⊻
        UInt64(item_id) * 0xbf58476d1ce4e5b9
    return _splitmix64(value)
end


"""Draw a categorical index from finite normalized or unnormalized masses."""
function categorical_draw(log_mass; rng::AbstractRNG)
    z = logsumexp(log_mass)
    isfinite(z) || throw(ArgumentError("proposal has no finite positive mass"))
    cumulative = cumsum(exp.(log_mass .- z))
    cumulative[end] = 1.0
    return something(findfirst(>(rand(rng)), cumulative))
end

"""Sample a virtual SIR population, preserving the pool mean-weight scale."""
function sir_pool_draws(log_mass, Ns, pool_count; rng::AbstractRNG)
    ids = [categorical_draw(log_mass; rng) for _ in 1:Ns]
    multiplicity = [count(==(i), ids) for i in eachindex(log_mass)]
    return (; selected_candidate_id=ids, multiplicity,
        active_candidate_id=findall(>(0), multiplicity),
        common_log_weight=logsumexp(log_mass)-log(pool_count))
end

"""Stop iid pool draws at Ns distinct IDs or 3Ns attempts, retaining repeats.

An undersized ID pool is retained with its original weights and multiplicities.
"""
function _sir_pool_draws_until_unique(log_weight, input_multiplicity, Ns; rng::AbstractRNG)
    Ns > 1 || throw(ArgumentError("SIR requires at least two target representatives"))
    pool_count = sum(input_multiplicity)
    _validate_multiplicity(log_weight, input_multiplicity, pool_count)
    all(m -> m isa Integer, input_multiplicity) || throw(ArgumentError("multiplicities must be integers"))
    K_pool = length(log_weight)
    ids = Int[]
    if K_pool < Ns
        return (; selected_candidate_id=ids, multiplicity=copy(input_multiplicity),
            active_candidate_id=collect(1:K_pool),
            candidate_output_log_weight=copy(log_weight),
            draw_count=0, stop_reason=:pool_under_target, K_pool,
            Nv_pool=pool_count, Nv_post=pool_count)
    end
    mass = log_weight .+ log.(input_multiplicity)
    log_total = logsumexp(mass)
    cumulative = cumsum(exp.(mass .- log_total))
    cumulative[end] = 1.0
    multiplicity = zeros(Int,K_pool)
    distinct = 0
    while distinct < Ns && length(ids) < 3Ns
        id = something(findfirst(>(rand(rng)),cumulative))
        push!(ids,id)
        distinct += iszero(multiplicity[id])
        multiplicity[id] += 1
    end
    return (; selected_candidate_id=ids, multiplicity,
        active_candidate_id=findall(>(0),multiplicity),
        candidate_output_log_weight=fill(log_total-log(pool_count),K_pool),
        draw_count=length(ids), stop_reason=distinct==Ns ? :target_reached : :draw_limit,
        K_pool, Nv_pool=pool_count, Nv_post=length(ids))
end

function _compensated_sum(values)
    total = 0.0
    correction = 0.0
    for value in values
        adjusted = value - correction
        updated = total + adjusted
        correction = (updated - total) - adjusted
        total = updated
    end
    return total
end

function sir_inclusion_probabilities(candidate_log_mass::AbstractVector{<:Real},
        number_draws::Integer)
    number_draws > 0 || throw(ArgumentError("SIR draw count must be positive"))
    all(x -> isfinite(x) || x == -Inf, candidate_log_mass) || throw(
        ArgumentError("candidate log masses must be finite or -Inf"),
    )
    active = findall(isfinite, candidate_log_mass)
    length(active) >= number_draws || throw(ArgumentError(
        "without-replacement SIR needs $number_draws positive candidate IDs; " *
        "found $(length(active))",
    ))
    shift = maximum(candidate_log_mass[active])
    shifted = Float64.(candidate_log_mass) .- shift
    order = sort(active; by=i -> shifted[i], rev=true, alg=MergeSort)
    suffix = zeros(Float64, length(order))
    suffix[end] = shifted[last(order)]
    for j in (length(order)-1):-1:1
        a, b = shifted[order[j]], suffix[j+1]
        hi, lo = max(a, b), min(a, b)
        suffix[j] = hi + log1p(exp(lo-hi))
    end
    log_total = first(suffix)
    log_threshold = -Inf # all positive candidates are retained
    if length(active) > number_draws
        h = 0
        while true
            log_threshold = suffix[h+1] - log(number_draws-h)
            shifted[order[h+1]] <= log_threshold && break
            h += 1
            h < number_draws || throw(ArgumentError("invalid SIR threshold"))
        end
    end
    log_inclusion = fill(-Inf, length(shifted))
    log_inclusion[active] .= min.(0.0, shifted[active] .- log_threshold)
    inclusion = exp.(log_inclusion)
    all(>(0), inclusion[active]) || throw(ArgumentError(
        "SIR inclusion probabilities underflowed; pool cannot be sampled in Float64",
    ))
    residual = _compensated_sum(inclusion) - number_draws
    abs(residual) <= 64eps(Float64)*number_draws || throw(ArgumentError(
        "SIR inclusion probabilities do not sum to the requested count",
    ))
    return (;
        probability=exp.(shifted .- log_total),
        inclusion_probability=inclusion,
        log_inclusion_probability=log_inclusion,
        certain_candidate_id=findall(==(1.0), inclusion),
        log_normalized_threshold=log_threshold-log_total,
        log_normalized_output_weight=[isfinite(x) ? max(x, log_threshold)-log_total :
            -Inf for x in shifted],
        # Add the common scale only after all relative arithmetic.
        log_total_mass=shift+log_total,
        inclusion_count_residual=residual,
    )
end

"""Select fixed-size, distinct IDs with the prescribed inclusion probabilities.
`order` must permute all candidates; `uniform` is one independent U[0,1) draw.
Separating this deterministic map also permits exact finite-pool validation.
"""
function sir_systematic_selection(inclusion::AbstractVector{<:Real},
        order::AbstractVector{<:Integer}, uniform::Real)
    all(x -> isfinite(x) && 0 <= x <= 1, inclusion) || throw(
        ArgumentError("inclusion probabilities must lie in [0,1]"),
    )
    sort(order) == collect(eachindex(inclusion)) || throw(ArgumentError(
        "systematic order must permute every candidate ID",
    ))
    isfinite(uniform) && 0 <= uniform < 1 || throw(ArgumentError(
        "systematic offset must lie in [0,1)",
    ))
    total = _compensated_sum(inclusion)
    n = round(Int, total)
    abs(total-n) <= 64eps(Float64)*max(1,n) || throw(ArgumentError(
        "systematic sampling requires an integer inclusion-probability sum",
    ))
    selected = Int[i for i in order if inclusion[i] == 1]
    light = Int[i for i in order if 0 < inclusion[i] < 1]
    remaining = n-length(selected)
    if remaining > 0
        cumulative = zeros(Float64, length(light))
        acc = correction = 0.0
        for (j,i) in enumerate(light)
            adjusted = inclusion[i]-correction
            updated = acc+adjusted
            correction = (updated-acc)-adjusted
            cumulative[j] = acc = updated
        end
        abs(acc-remaining) <= 64eps(Float64)*max(1,n) || throw(
            ArgumentError("systematic cumulative endpoint disagrees with count"),
        )
        cumulative[end] = Float64(remaining)
        index = 1
        for j in 0:(remaining-1)
            point = uniform+j
            while index <= length(light) && cumulative[index] <= point
                index += 1
            end
            index <= length(light) || throw(ArgumentError(
                "systematic sample missed its cumulative interval",
            ))
            push!(selected, light[index])
        end
    end
    length(selected) == n && allunique(selected) || throw(ArgumentError(
        "systematic sample did not select exactly n distinct candidate IDs",
    ))
    return selected
end

"""
Without-replacement SIR using capped inclusion probabilities and HT weights.
Two explicit RNGs independently choose a permutation and systematic offset.
Output log weights preserve `sum(pool masses)/pool_virtual_count` as the
ensemble mean weight. The stored Born proposal and physical norm are unrelated
to the conditional inclusion probability and must be retained by the caller.
"""
function sir_pool_draws_without_replacement(candidate_log_mass::AbstractVector{<:Real},
        number_draws::Integer, pool_virtual_count::Integer;
        order_rng::AbstractRNG, offset_rng::AbstractRNG)
    pool_virtual_count >= length(candidate_log_mass) || throw(ArgumentError(
        "virtual pool count cannot be below its representative count",
    ))
    design = sir_inclusion_probabilities(candidate_log_mass, number_draws)
    order = randperm(order_rng, length(candidate_log_mass))
    uniform = rand(offset_rng)
    selected = sir_systematic_selection(design.inclusion_probability, order, uniform)
    multiplicity = zeros(Int, length(candidate_log_mass))
    multiplicity[selected] .= 1
    active = findall(>(0), multiplicity)
    normalized_output = exp.(design.log_normalized_output_weight[active])
    residual = _compensated_sum(normalized_output)-1
    abs(residual) <= 128eps(Float64)*max(1,number_draws) || throw(
        ArgumentError("without-replacement SIR did not preserve total mass"),
    )
    log_scale = design.log_total_mass + log(number_draws)-log(pool_virtual_count)
    certain = design.certain_candidate_id
    second_moment = _compensated_sum(exp.(
        2 .* design.log_normalized_output_weight[certain])) +
        (number_draws-length(certain))*exp(2design.log_normalized_threshold)
    return (;
        selected_candidate_id=selected,
        multiplicity,
        active_candidate_id=active,
        candidate_normalized_mass=design.probability,
        candidate_output_log_weight=log_scale .+ design.log_normalized_output_weight,
        candidate_inclusion_probability=design.inclusion_probability,
        candidate_log_inclusion_probability=design.log_inclusion_probability,
        certain_candidate_id=design.certain_candidate_id,
        log_normalized_threshold=design.log_normalized_threshold,
        order,
        uniform,
        inclusion_count_residual=design.inclusion_count_residual,
        normalized_mass_residual=residual,
        eta_expected=inv(number_draws*second_moment),
    )
end
