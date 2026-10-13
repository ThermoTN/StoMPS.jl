function state_bond_dimensions(state)
    return Int[
        dim(domain(state[site].A)[end])
        for site in 1:max(length(state) - 1, 0)
    ]
end

function state_max_bond_dimension(state)
    dimensions = state_bond_dimensions(state)
    return isempty(dimensions) ? 1 : maximum(dimensions)
end

function tdvp_second_half_diagnostics(info, state)
    number_bonds = max(length(state) - 1, 0)
    entropy = zeros(Float64, number_bonds)
    truncation = zeros(Float64, number_bonds)
    for bond in 1:number_bonds
        bond_info = info[2].backward[bond].Bond
        entropy[bond] = bond_info.SE
        truncation[bond] = bond_info.TrunErr^2
    end
    dimensions = state_bond_dimensions(state)
    return (
        entropy,
        truncation,
        bond_dimension=dimensions,
        maximum_entropy=isempty(entropy) ? 0.0 : maximum(entropy),
        maximum_truncation=isempty(truncation) ? 0.0 : maximum(truncation),
        maximum_bond_dimension=isempty(dimensions) ? 1 : maximum(dimensions),
    )
end
