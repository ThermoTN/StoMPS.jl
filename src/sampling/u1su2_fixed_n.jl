"""
    u1su2_born_space()

Construct the canonical local Hubbard space for the
U(1)×SU(2) Born backend. The U(1) charge is `q=n-1`, corresponding to
`C_z-1` in Fig. 1(e) of [arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1).
The public measurement
order is empty, single, double occupancy, with sectors
`(-1,0)`, `(0,1/2)`, `(1,0)` and carrier dimensions `1,2,1`.
"""
function u1su2_born_space()
    return Rep[U₁ × SU₂]((-1, 0) => 1, (0, 1 // 2) => 1, (1, 0) => 1)
end

function _u1su2_sector_list()
    return [
        Irrep[U₁ × SU₂](-1, 0),
        Irrep[U₁ × SU₂](0, 1 // 2),
        Irrep[U₁ × SU₂](1, 0),
    ]
end

_u1su2_occupancy(sector) = Int(sector[1].charge) + 1

function _validate_u1su2_space(physical_space::GradedSpace)
    physical_space == u1su2_born_space() || throw(ArgumentError(
        "U(1)×SU(2) Born sampling requires canonical u1su2_born_space()",
    ))
    sector_list = _u1su2_sector_list()
    all(sector -> dim(physical_space, sector) == 1, sector_list) || throw(
        ArgumentError("U(1)×SU(2) Born sampling requires unit multiplicity"),
    )
    marker_dimensions = [
        dim(_born_sector_space(physical_space, sector)) for sector in sector_list
    ]
    marker_dimensions == [1, 2, 1] || throw(ArgumentError(
        "U(1)×SU(2) Born markers must have dimensions [1,2,1]",
    ))
    return (; sector_list, marker_dimensions)
end

function _u1su2_projected_tensor(tensor, physical_space, sector)
    marker_space = _born_sector_space(physical_space, sector)
    embedding = isometry(physical_space, marker_space)
    local_tensor = tensor.A
    @tensor projected[left physical; marker right] :=
        local_tensor[left physical; input right] * embedding[input; marker]
    return MPSTensor(projected)
end

# Identity-channel transfers for rank-4 purification tensors, using public
# FiniteMPS tensor wrappers and TensorKit contractions.
function _u1su2_pushleft(environment::LocalRightTensor{2}, tensor::MPSTensor{4})
    ket = tensor.A
    bra = adjoint(tensor).A
    @tensor next[left; left_prime] :=
        (ket[left physical marker right] * environment.A[right; right_prime]) *
        bra[marker right_prime left_prime physical]
    return LocalRightTensor(next, environment.tag)
end

function _u1su2_pushright(environment::LocalLeftTensor{2}, tensor::MPSTensor{4})
    ket = tensor.A
    bra = adjoint(tensor).A
    @tensor next[right_prime; right] :=
        (environment.A[left_prime; left] *
         bra[marker right_prime left_prime physical]) *
        ket[left physical marker right]
    return LocalLeftTensor(next, environment.tag)
end

function _u1su2_add_environment!(dictionary, charge::Int, contribution)
    if haskey(dictionary, charge)
        current = dictionary[charge]
        dictionary[charge] = LocalRightTensor(
            current.A + contribution.A,
            current.tag,
        )
    else
        dictionary[charge] = contribution
    end
    return nothing
end

function _u1su2_projected_tensors(thermal_mpo, physical_space, sector_list)
    return [
        [_u1su2_projected_tensor(thermal_mpo[site], physical_space, sector)
         for sector in sector_list]
        for site in eachindex(thermal_mpo)
    ]
end

function _u1su2_right_charge_environments(projected, boundary, charges)
    L = length(projected)
    right = [Dict{Int,LocalRightTensor{2}}() for _ in 1:(L + 1)]
    right[L + 1][0] = boundary
    for site in L:-1:1
        for (suffix_charge, suffix_environment) in right[site + 1]
            for sector_index in eachindex(charges)
                contribution = _u1su2_pushleft(
                    suffix_environment,
                    projected[site][sector_index],
                )
                _u1su2_add_environment!(
                    right[site],
                    charges[sector_index] + suffix_charge,
                    contribution,
                )
            end
        end
    end
    return right
end

function _u1su2_fixed_n_select(
    thermal_mpo::MPO{L},
    physical_space::GradedSpace,
    target_N::Integer;
    select_sector,
    trunc=notrunc(),
    probability_tolerance::Real=1e-10,
) where {L}
    probability_tolerance > 0 || throw(ArgumentError(
        "probability_tolerance must be positive",
    ))
    0 <= target_N <= 2L || throw(ArgumentError(
        "target_N must satisfy 0 <= target_N <= 2L",
    ))
    validated = _validate_u1su2_space(physical_space)
    sector_list = validated.sector_list
    charges = Int[sector[1].charge for sector in sector_list]
    occupancies_in_order = charges .+ 1
    for site in eachindex(thermal_mpo)
        tensor = thermal_mpo[site]
        tensor isa MPSTensor{4} || throw(ArgumentError(
            "fixed-N Born sampling requires rank-4 tensors",
        ))
        codomain(tensor)[2] == physical_space || throw(ArgumentError(
            "thermal MPO output physical-space mismatch at site $site",
        ))
        domain(tensor)[1] == physical_space || throw(ArgumentError(
            "thermal MPO input physical-space mismatch at site $site",
        ))
    end
    isapprox(norm(thermal_mpo), 1; atol=probability_tolerance, rtol=0) || throw(
        ArgumentError("fixed-N Born sampling requires a normalized purification"),
    )

    projected = _u1su2_projected_tensors(thermal_mpo, physical_space, sector_list)
    environment = Environment(thermal_mpo', thermal_mpo; disk=false)
    right = _u1su2_right_charge_environments(projected, environment.Er[L], charges)
    target_charge = Int(target_N) - L
    haskey(right[1], target_charge) || throw(ArgumentError(
        "target particle-number sector has zero support",
    ))
    target_sector_weight = abs2(FiniteMPS.coef(thermal_mpo)) *
        real(environment.El[1] * right[1][target_charge])
    isfinite(target_sector_weight) && target_sector_weight > 0 || throw(
        ArgumentError("target particle-number sector has zero or invalid weight"),
    )

    left = environment.El[1]
    configuration = Int[]
    occupancies = Int[]
    conditional_probabilities = Float64[]
    probability_residuals = Float64[]
    prefix_charge = 0
    for site in 1:L
        raw_weights = zeros(Float64, length(sector_list))
        for sector_index in eachindex(sector_list)
            remaining = target_charge - prefix_charge - charges[sector_index]
            haskey(right[site + 1], remaining) || continue
            candidate_left = _u1su2_pushright(
                left,
                projected[site][sector_index],
            )
            value = real(candidate_left * right[site + 1][remaining])
            value >= -probability_tolerance || error(
                "fixed-N Born sampling produced a negative branch weight at site $site",
            )
            raw_weights[sector_index] = max(value, 0.0)
        end
        total = sum(raw_weights)
        isfinite(total) && total > 0 || error(
            "fixed-N Born probability vanished at site $site",
        )
        probabilities = raw_weights ./ total
        residual = abs(sum(probabilities) - 1)
        residual <= probability_tolerance || error(
            "fixed-N Born probability normalization failed at site $site",
        )
        selected = select_sector(probabilities, site)
        selected isa Integer && selected in eachindex(sector_list) || throw(
            ArgumentError("Born selector returned an invalid sector at site $site"),
        )
        probabilities[selected] > 0 || throw(ArgumentError(
            "selected fixed-N configuration has zero probability at site $site",
        ))
        push!(configuration, selected)
        push!(occupancies, occupancies_in_order[selected])
        push!(conditional_probabilities, probabilities[selected])
        push!(probability_residuals, residual)
        prefix_charge += charges[selected]
        left = _u1su2_pushright(left, projected[site][selected])
    end
    prefix_charge == target_charge || error("fixed-N sampler ended in the wrong charge")

    selected_tensors = AbstractMPSTensor[
        projected[site][configuration[site]] for site in 1:L
    ]
    state = MPS(selected_tensors, [1, L], FiniteMPS.coef(thermal_mpo))
    canonicalize!(state, L; trunc=notrunc())
    canonicalize!(state, 1; trunc=notrunc())
    projection_norm = norm(state)
    isfinite(projection_norm) && projection_norm > 0 || error(
        "fixed-N projection produced a zero or invalid state",
    )
    unconditioned_born_probability = abs2(projection_norm)
    normalize!(state)
    canonicalize!(state, 1; trunc=trunc)
    compression_norm = norm(state)
    isfinite(compression_norm) && compression_norm > 0 || error(
        "fixed-N projection compression produced a zero or invalid state",
    )
    normalize!(state)

    born_probability = unconditioned_born_probability / target_sector_weight
    log_probability = sum(log, conditional_probabilities)
    isapprox(exp(log_probability), born_probability;
        atol=10probability_tolerance, rtol=10probability_tolerance) || error(
        "sequential and direct fixed-N Born probabilities disagree",
    )
    return (
        state,
        configuration,
        occupancies,
        conditional_probabilities,
        log_probability,
        born_probability,
        unconditioned_born_probability,
        target_sector_weight,
        target_N=Int(target_N),
        target_charge,
        sector_list,
        marker_dimensions=validated.marker_dimensions,
        compression_norm,
        maximum_probability_residual=maximum(probability_residuals),
    )
end

"""
    u1su2_fixed_n_purification_born_sample(thermal_mpo, target_N; rng, ...)

Draw a normalized fixed-particle-number U(1)×SU(2) Born sample. The
returned `log_probability` is the conditional law `log p_beta0,Ne(alpha)`,
as defined in App. D, Eq. (24) of
[arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1).
"""
function u1su2_fixed_n_purification_born_sample(
    thermal_mpo::MPO{L},
    target_N::Integer;
    rng::AbstractRNG,
    trunc=notrunc(),
    probability_tolerance::Real=1e-10,
) where {L}
    physical_space = codomain(thermal_mpo[1])[2]
    return _u1su2_fixed_n_select(
        thermal_mpo,
        physical_space,
        target_N;
        select_sector=(probabilities, _site) ->
            _draw_born_probability_index(probabilities, rng),
        trunc,
        probability_tolerance,
    )
end

"""
    u1su2_fixed_n_purification_born_project(
        thermal_mpo, target_N, configuration; ...)

Deterministically project the normalized purification onto a complete local
occupancy-sector configuration in the fixed-`target_N` conditional law.
"""
function u1su2_fixed_n_purification_born_project(
    thermal_mpo::MPO{L},
    target_N::Integer,
    configuration::AbstractVector{<:Integer};
    trunc=notrunc(),
    probability_tolerance::Real=1e-10,
) where {L}
    length(configuration) == L || throw(DimensionMismatch(
        "configuration length $(length(configuration)) does not equal $L",
    ))
    all(index -> index in 1:3, configuration) || throw(ArgumentError(
        "U(1)×SU(2) configuration contains an invalid sector index",
    ))
    sum(index -> index - 1, configuration) == target_N || throw(ArgumentError(
        "configuration particle number does not equal target_N",
    ))
    physical_space = codomain(thermal_mpo[1])[2]
    return _u1su2_fixed_n_select(
        thermal_mpo,
        physical_space,
        target_N;
        select_sector=(_probabilities, site) -> configuration[site],
        trunc,
        probability_tolerance,
    )
end
