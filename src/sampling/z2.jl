"""
    z2_born_space(charges::Integer...)

Construct a dictionary-backed, unit-multiplicity Z2 graded space for the
selected charges.  Charges must be `0` or `1` and may not be repeated.

This constructor keeps the physical and auxiliary spaces used by the first
StoMPS Born backend on one TensorKit storage representation.  For example,
`z2_born_space(0, 1)` is the supported local physical space and
`z2_born_space(1)` is the odd auxiliary space used by the transverse-field
Ising interaction.
"""
function z2_born_space(charges::Integer...)
    isempty(charges) && throw(ArgumentError("at least one Z2 charge is required"))
    all(charge -> charge == 0 || charge == 1, charges) || throw(ArgumentError(
        "Z2 charges must be 0 or 1",
    ))
    allunique(charges) || throw(ArgumentError("Z2 charges may not be repeated"))

    sector_pairs = map(charges) do charge
        Irrep[ℤ₂](charge) => 1
    end
    dimensions = FiniteMPS.TensorKit.SortedVectorDict(sector_pairs...)
    sector_type = typeof(Irrep[ℤ₂](0))
    return GradedSpace{sector_type,typeof(dimensions)}(dimensions, false)
end

"""
    z2_purification_born_sample(thermal_mpo, physical_space; rng, trunc=notrunc())

Sequentially measure the `physical_in` legs of a normalized rank-4 thermal MPO
in the local Z2 charge basis.  Interpreting

```text
thermal_mpo[left, physical_out; physical_in, right]
```

as a purification, this samples a configuration `a` with the exact Born
probability

```text
p(a) = norm(thermal_mpo * |a>)^2
```

and returns the normalized physical snapshot

```text
thermal_mpo * |a> / sqrt(p(a)).
```

The returned state keeps a one-dimensional graded marker in the third local
leg, so it remains a rank-4 `MPS` compatible with the existing StoMPS TDVP
workflow.  Its virtual bonds are inherited from the thermal MPO.  Local Z2
sector multiplicities must be one.
"""
function z2_purification_born_sample(
    thermal_mpo::MPO{L},
    physical_space::GradedSpace;
    rng::AbstractRNG,
    trunc=notrunc(),
    probability_tolerance::Float64=1e-10,
) where {L}
    sector_list = collect(sectors(physical_space))
    expected_sectors = Set((Irrep[ℤ₂](0), Irrep[ℤ₂](1)))
    Set(sector_list) == expected_sectors || throw(ArgumentError(
        "Z2 Born sampling requires exactly the charge-0 and charge-1 sectors",
    ))
    all(sector -> dim(physical_space, sector) == 1, sector_list) || throw(
        ArgumentError("Z2 Born sampling currently requires unit sector multiplicities"),
    )
    all(sector -> dim(_born_sector_space(physical_space, sector)) == 1, sector_list) ||
        throw(ArgumentError("Z2 Born markers must be one dimensional"))
    return _purification_sector_born_sample(
        thermal_mpo,
        physical_space,
        sector_list;
        rng,
        trunc,
        probability_tolerance,
    )
end
