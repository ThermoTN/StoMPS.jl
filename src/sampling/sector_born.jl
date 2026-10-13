function _purification_sector_born_select(
    thermal_mpo::MPO{L},
    physical_space::GradedSpace,
    sector_list::Vector;
    select_sector::Selector,
    trunc=notrunc(),
    probability_tolerance::Float64=1e-10,
) where {L,Selector}
    L > 0 || throw(ArgumentError("Born sampling requires a nonempty thermal MPO"))
    isapprox(norm(thermal_mpo), 1; atol=probability_tolerance, rtol=0) || throw(
        ArgumentError("Born sampling requires a normalized thermal MPO"),
    )
    probability_tolerance > 0 || throw(ArgumentError(
        "probability_tolerance must be positive",
    ))
    isempty(sector_list) && throw(ArgumentError("physical space has no sectors"))
    allunique(sector_list) || throw(ArgumentError("Born sectors may not be repeated"))
    all(sector -> sector in sectors(physical_space), sector_list) || throw(
        ArgumentError("Born sector is absent from the physical space"),
    )

    marker_spaces = map(sector_list) do sector
        _born_sector_space(physical_space, sector)
    end
    tensors = AbstractMPSTensor[deepcopy(tensor) for tensor in thermal_mpo]
    snapshot = MPS(
        tensors,
        copy(Center(thermal_mpo)),
        FiniteMPS.coef(thermal_mpo),
    )
    canonicalize!(snapshot, 1; trunc=notrunc())
    normalize!(snapshot)

    configuration = Vector{Int}(undef, L)
    conditional_probabilities = Vector{Float64}(undef, L)
    probability_residuals = Vector{Float64}(undef, L)
    log_probability = 0.0

    for site in 1:L
        tensor = snapshot[site]
        tensor isa MPSTensor{4} || throw(ArgumentError(
            "Born sampling requires a rank-4 tensor at site $site",
        ))
        codomain(tensor)[2] == physical_space || throw(ArgumentError(
            "Born sampling output physical-space mismatch at site $site",
        ))
        domain(tensor)[1] == physical_space || throw(ArgumentError(
            "Born sampling input physical-space mismatch at site $site",
        ))

        candidates = map(marker_spaces) do marker_space
            embedding = isometry(physical_space, marker_space)
            local_tensor = tensor.A
            @tensor projected[left physical; marker right] :=
                local_tensor[left physical; input right] *
                embedding[input; marker]
            MPSTensor(projected)
        end
        raw_probabilities = abs2.(norm.(candidates))
        probability_sum = sum(raw_probabilities)
        isfinite(probability_sum) && probability_sum > 0 || error(
            "Born sampling probability vanished at site $site",
        )
        probability_residuals[site] = abs(probability_sum - 1)
        probability_residuals[site] <= probability_tolerance || error(
            "Born probability normalization failed at site $site: " *
            "sum=$probability_sum, residual=$(probability_residuals[site])",
        )
        probabilities = raw_probabilities ./ probability_sum

        selected_index = select_sector(probabilities, site)
        selected_index isa Integer && selected_index in eachindex(probabilities) || throw(
            ArgumentError("Born selector returned invalid sector index $selected_index at site $site"),
        )
        selected_probability = probabilities[selected_index]
        selected_probability > 0 || error(
            "Born sampling selected a zero-probability sector at site $site",
        )
        configuration[site] = selected_index
        conditional_probabilities[site] = selected_probability
        log_probability += log(selected_probability)

        selected_norm = norm(candidates[selected_index])
        snapshot[site] = candidates[selected_index] / selected_norm
        if site < L
            canonicalize!(snapshot, site + 1)
            conditional_norm = norm(snapshot)
            isfinite(conditional_norm) && conditional_norm > 0 || error(
                "Born snapshot became zero or non-finite after site $site",
            )
            normalize!(snapshot)
        end
    end

    canonicalize!(snapshot, 1; trunc=trunc)
    compression_norm = norm(snapshot)
    isfinite(compression_norm) && compression_norm > 0 || error(
        "Born snapshot compression produced a zero or non-finite norm",
    )
    normalize!(snapshot)

    return (
        state=snapshot,
        configuration,
        conditional_probabilities,
        log_probability,
        compression_norm,
        maximum_probability_residual=maximum(probability_residuals),
        sector_list,
    )
end

function _purification_sector_born_sample(
    thermal_mpo::MPO{L},
    physical_space::GradedSpace,
    sector_list::Vector;
    rng::AbstractRNG,
    trunc=notrunc(),
    probability_tolerance::Float64=1e-10,
) where {L}
    return _purification_sector_born_select(
        thermal_mpo,
        physical_space,
        sector_list;
        select_sector=(probabilities, _site) ->
            _draw_born_probability_index(probabilities, rng),
        trunc,
        probability_tolerance,
    )
end

function _draw_born_probability_index(probabilities, rng::AbstractRNG)
    threshold = rand(rng)
    cumulative_probability = 0.0
    selected_index = lastindex(probabilities)
    for index in eachindex(probabilities)
        cumulative_probability += probabilities[index]
        if threshold < cumulative_probability
            selected_index = index
            break
        end
    end
    return selected_index
end

function _born_sector_space(
    physical_space::GradedSpace{Sector,Storage},
    sector::Sector,
) where {Sector,Storage}
    dimensions = FiniteMPS.TensorKit.SortedVectorDict(
        sector => dim(physical_space, sector),
    )
    typeof(dimensions) == Storage || throw(ArgumentError(
        "Born marker storage $(typeof(dimensions)) does not match physical-space " *
        "storage $Storage",
    ))
    return GradedSpace{Sector,Storage}(dimensions, false)
end
