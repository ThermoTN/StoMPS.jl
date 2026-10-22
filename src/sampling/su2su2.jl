"""
    su2su2_born_space()

Construct the unit-multiplicity local space for the half-filled Hubbard SU(2)×SU(2)
backend.  Its two factors are ordered `(charge SU(2), spin SU(2))` and its
admitted irreps are `(1/2,0)` and `(0,1/2)`.
"""
function su2su2_born_space()
    return Rep[SU₂ × SU₂]((1 // 2, 0) => 1, (0, 1 // 2) => 1)
end

"""
    su2su2_purification_born_sample(thermal_mpo, physical_space;
        rng, trunc=notrunc())

Sequentially measure the `physical_in` legs of a normalized SU(2)charge × SU(2)spin
thermal MPO with the two local irreducible-sector projectors

```text
P_eta  : (C=1/2, S=0)     empty/double-occupancy charge doublet
P_spin : (C=0,   S=1/2)   single-occupancy spin doublet.
```

This symmetry-preserving Born measurement samples the local multiplet label,
as in Sec. II B of [arXiv:2609.35754v1](https://arxiv.org/html/2609.35754v1).
The selected two-dimensional multiplet is retained as the third (`marker`) leg
of each rank-4 sample tensor. Consequently every normalized Born sample
preserves charge SU(2) and spin SU(2), and the two projectors satisfy
`P_eta + P_spin = I` locally.

The physical space must be exactly
`Rep[SU2 x SU2]((1/2,0)=>1, (0,1/2)=>1)`, with factors ordered as
`(charge, spin)`.  The return value extends the common Born diagnostics with
`marker_dimensions`, which is `[2,2]` in `sector_list` order.
"""
function su2su2_purification_born_sample(
    thermal_mpo::MPO{L},
    physical_space::GradedSpace;
    rng::AbstractRNG,
    trunc=notrunc(),
    probability_tolerance::Float64=1e-10,
) where {L}
    sector_list = collect(sectors(physical_space))
    canonical_space = su2su2_born_space()
    expected_sectors = Set((
        Irrep[SU₂ × SU₂](1 // 2, 0),
        Irrep[SU₂ × SU₂](0, 1 // 2),
    ))
    physical_space == canonical_space || throw(ArgumentError(
        "SU2SU2 Born sampling requires the canonical su2su2_born_space()",
    ))
    Set(sector_list) == expected_sectors || throw(ArgumentError(
        "SU2SU2 Born sampling requires exactly the (1/2,0) and (0,1/2) sectors",
    ))
    all(sector -> dim(physical_space, sector) == 1, sector_list) || throw(
        ArgumentError("SU2SU2 Born sampling requires unit sector multiplicities"),
    )
    marker_dimensions = map(sector_list) do sector
        dim(_born_sector_space(physical_space, sector))
    end
    marker_dimensions == [2, 2] || throw(ArgumentError(
        "SU2SU2 Born markers must both be two dimensional; got $marker_dimensions",
    ))

    result = _purification_sector_born_sample(
        thermal_mpo,
        physical_space,
        sector_list;
        rng,
        trunc,
        probability_tolerance,
    )
    return merge(result, (; marker_dimensions))
end

"""
    su2su2_purification_born_project(thermal_mpo, physical_space, configuration;
                                  trunc=notrunc())

Deterministically project a normalized SU(2)charge × SU(2)spin thermal
purification onto a full charge/spin multiplet configuration. This is the
fixed-selector counterpart of
[`su2su2_purification_born_sample`](@ref), and shares its sequential contractions,
auxiliary-leg convention, normalization checks, and final compression policy.
The selected configuration must have positive Born probability.
"""
function su2su2_purification_born_project(
    thermal_mpo::MPO{L},
    physical_space::GradedSpace,
    configuration::AbstractVector{<:Integer};
    trunc=notrunc(),
    probability_tolerance::Float64=1e-10,
) where {L}
    length(configuration) == L || throw(DimensionMismatch(
        "SU2SU2 configuration requires length $L; got $(length(configuration))",
    ))
    sector_list = collect(sectors(physical_space))
    _validate_su2su2_resampling_space(physical_space)
    all(index -> index in eachindex(sector_list), configuration) || throw(
        ArgumentError("SU2SU2 configuration contains an invalid sector index"),
    )
    marker_dimensions = [
        dim(_born_sector_space(physical_space, sector)) for sector in sector_list
    ]
    result = _purification_sector_born_select(
        thermal_mpo,
        physical_space,
        sector_list;
        select_sector=(_probabilities, site) -> configuration[site],
        trunc,
        probability_tolerance,
    )
    return merge(result, (; marker_dimensions))
end
