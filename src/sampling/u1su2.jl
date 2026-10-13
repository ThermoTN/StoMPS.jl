function _u1su2_unconditioned_select(thermal_mpo::MPO{L}; select_sector,
        trunc=notrunc(), probability_tolerance::Real=1e-10) where {L}
    L > 0 || throw(ArgumentError("Born sampling requires a nonempty purification"))
    isfinite(probability_tolerance) && probability_tolerance > 0 ||
        throw(ArgumentError("probability_tolerance must be finite and positive"))
    physical_space = codomain(thermal_mpo[1])[2]
    validated = _validate_u1su2_space(physical_space)
    result = _purification_sector_born_select(thermal_mpo, physical_space,
        validated.sector_list; select_sector, trunc,
        probability_tolerance=Float64(probability_tolerance))
    occupancies = result.configuration .- 1
    particle_number = sum(occupancies)
    return merge(result, (; occupancies, particle_number,
        total_charge=particle_number - L,
        born_probability=exp(result.log_probability),
        marker_dimensions=validated.marker_dimensions))
end

"""
    u1su2_purification_born_sample(thermal_mpo; rng, trunc=notrunc(),
        probability_tolerance=1e-10)

Sample the **unconditioned** Born law of a normalized U(1)×SU(2) purification.
Measure complete empty/single/double projectors in `u1su2_born_space()` order,
retaining carrier dimensions `[1,2,1]` as rank-4 marker legs. No magnetic spin
quantum number is resolved. Each snapshot has a definite physical particle
number, but different snapshots may have different particle numbers.

Return `state`, `configuration`, `occupancies`, `particle_number`,
`total_charge=particle_number-length(thermal_mpo)`, `conditional_probabilities`,
`log_probability`, `born_probability`, `sector_list`, `marker_dimensions`,
`compression_norm`, and `maximum_probability_residual`. Probabilities describe
the uncompressed projection. The returned state is normalized after compression;
the compression norm is a numerical diagnostic, not an importance weight.
The exact snapshot convention is `trunc=notrunc()`.
"""
function u1su2_purification_born_sample(thermal_mpo::MPO;
        rng::AbstractRNG, trunc=notrunc(), probability_tolerance::Real=1e-10)
    return _u1su2_unconditioned_select(thermal_mpo;
        select_sector=(probabilities, _site) ->
            _draw_born_probability_index(probabilities, rng),
        trunc, probability_tolerance)
end

"""
    u1su2_purification_born_project(thermal_mpo, configuration;
        trunc=notrunc(), probability_tolerance=1e-10)

Deterministically project onto a complete empty/single/double configuration
(indices `1:3`) using the unconditioned Born law. Return the same fields as
[`u1su2_purification_born_sample`](@ref). A zero-probability configuration is
rejected; this operation neither conditions on nor requests a target N.
"""
function u1su2_purification_born_project(thermal_mpo::MPO{L},
        configuration::AbstractVector{<:Integer}; trunc=notrunc(),
        probability_tolerance::Real=1e-10) where {L}
    length(configuration) == L || throw(DimensionMismatch("configuration length differs from purification"))
    all(in(1:3), configuration) || throw(ArgumentError("configuration indices must lie in 1:3"))
    return _u1su2_unconditioned_select(thermal_mpo;
        select_sector=(_probabilities, site) -> configuration[site],
        trunc, probability_tolerance)
end
