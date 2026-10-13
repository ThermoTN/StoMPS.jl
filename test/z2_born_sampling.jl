using Random
using StoMPS.FiniteMPS
using StoMPS.FiniteMPS.TensorKit

mutable struct ThresholdRNG <: AbstractRNG
    values::Vector{Float64}
    position::Int
end

ThresholdRNG(values) = ThresholdRNG(collect(values), 1)

function z2_charge(sector)
    sector == Irrep[ℤ₂](0) && return 0
    sector == Irrep[ℤ₂](1) && return 1
    throw(ArgumentError("not a Z2 charge-0 or charge-1 sector"))
end

function Random.rand(rng::ThresholdRNG)
    rng.position <= length(rng.values) || error("ThresholdRNG is exhausted")
    value = rng.values[rng.position]
    rng.position += 1
    return value
end

function reference_thermal_mpo(physical_space)
    boundary_space = z2_born_space(0)
    virtual_space = z2_born_space(0, 1)
    charge_zero = Irrep[ℤ₂](0)
    charge_one = Irrep[ℤ₂](1)

    first_tensor = zeros(
        Float64,
        boundary_space ⊗ physical_space,
        physical_space ⊗ virtual_space,
    )
    block(first_tensor, charge_zero) .= [1.0 0.2]
    block(first_tensor, charge_one) .= [0.7 -0.3]

    second_tensor = zeros(
        Float64,
        virtual_space ⊗ physical_space,
        physical_space ⊗ boundary_space,
    )
    block(second_tensor, charge_zero) .= [1.1; 0.4]
    block(second_tensor, charge_one) .= [0.6; -0.25]

    thermal_mpo = MPO([MPSTensor(first_tensor), MPSTensor(second_tensor)])
    canonicalize!(thermal_mpo, length(thermal_mpo); trunc=notrunc())
    canonicalize!(thermal_mpo, 1; trunc=notrunc())
    normalize!(thermal_mpo)
    return thermal_mpo
end

function direct_column(thermal_mpo, physical_space, configuration)
    sector_list = collect(sectors(physical_space))
    tensors = AbstractMPSTensor[]
    for site in eachindex(configuration)
        sector = sector_list[configuration[site]]
        marker_space = z2_born_space(z2_charge(sector))
        embedding = isometry(physical_space, marker_space)
        local_tensor = thermal_mpo[site].A
        @tensor projected[left physical; marker right] :=
            local_tensor[left physical; input right] *
            embedding[input; marker]
        push!(tensors, MPSTensor(projected))
    end

    column = MPS(tensors, [1, length(tensors)], FiniteMPS.coef(thermal_mpo))
    canonicalize!(column, length(column); trunc=notrunc())
    canonicalize!(column, 1; trunc=notrunc())
    probability = abs2(norm(column))
    probability > 0 && normalize!(column)
    return (; state=column, probability)
end

function basis_state(physical_space, output_configuration, input_configuration)
    sector_list = collect(sectors(physical_space))
    tensors = AbstractMPSTensor[]
    left_charge = 0

    for site in eachindex(output_configuration)
        output_charge = z2_charge(sector_list[output_configuration[site]])
        input_charge = z2_charge(sector_list[input_configuration[site]])
        right_charge = mod(left_charge + output_charge + input_charge, 2)

        left_space = z2_born_space(left_charge)
        input_marker = z2_born_space(input_charge)
        right_space = z2_born_space(right_charge)
        tensor = ones(
            Float64,
            left_space ⊗ physical_space,
            input_marker ⊗ right_space,
        )
        push!(tensors, MPSTensor(tensor))
        left_charge = right_charge
    end

    left_charge == 0 || throw(ArgumentError(
        "output and input configurations must have the same total Z2 charge",
    ))
    state = MPS(tensors)
    canonicalize!(state, length(state); trunc=notrunc())
    canonicalize!(state, 1; trunc=notrunc())
    normalize!(state)
    return state
end

function configuration_thresholds(configuration, configurations, probabilities)
    thresholds = Float64[]
    for site in eachindex(configuration)
        prefix_matches(candidate) = all(
            candidate[index] == configuration[index] for index in 1:(site - 1)
        )
        denominator = sum(
            probabilities[index] for index in eachindex(configurations)
            if prefix_matches(configurations[index])
        )
        local_probabilities = map(1:2) do sector_index
            sum(
                probabilities[index] for index in eachindex(configurations)
                if prefix_matches(configurations[index]) &&
                   configurations[index][site] == sector_index
            ) / denominator
        end
        selected = configuration[site]
        push!(
            thresholds,
            sum(local_probabilities[1:(selected - 1)]; init=0.0) +
            local_probabilities[selected] / 2,
        )
    end
    return thresholds
end

function dense_tfim_hamiltonian(configurations; coupling=0.8, field=0.35)
    dimension = length(configurations)
    hamiltonian = zeros(Float64, dimension, dimension)
    configuration_index = Dict(config => index for (index, config) in enumerate(configurations))

    for (index, configuration) in enumerate(configurations)
        sigma_x = map(charge_index -> charge_index == 1 ? 1.0 : -1.0, configuration)
        hamiltonian[index, index] -= field * sum(sigma_x)
        flipped = (3 - configuration[1], 3 - configuration[2])
        hamiltonian[configuration_index[flipped], index] -= coupling
    end
    return hamiltonian
end

@testset "Z2 Born spaces" begin
    physical_space = z2_born_space(0, 1)
    @test collect(sectors(physical_space)) == [Irrep[ℤ₂](0), Irrep[ℤ₂](1)]
    @test all(sector -> dim(physical_space, sector) == 1, sectors(physical_space))
end

@testset "Z₂ purification Born sampling" begin
    physical_space = z2_born_space(0, 1)
    number_sites = 3

    unnormalized_mpo = identityMPO(number_sites, physical_space)
    @test_throws ArgumentError z2_purification_born_sample(
        unnormalized_mpo,
        physical_space;
        rng=MersenneTwister(1),
    )

    wrong_symmetry_mpo = identityMPO(number_sites, SU2Spin.pspace)
    normalize!(wrong_symmetry_mpo)
    @test_throws ArgumentError z2_purification_born_sample(
        wrong_symmetry_mpo,
        SU2Spin.pspace;
        rng=MersenneTwister(1),
    )

    thermal_mpo = deepcopy(unnormalized_mpo)
    normalize!(thermal_mpo)
    first_sample = z2_purification_born_sample(
        thermal_mpo,
        physical_space;
        rng=MersenneTwister(17),
        trunc=notrunc(),
    )
    repeated_sample = z2_purification_born_sample(
        thermal_mpo,
        physical_space;
        rng=MersenneTwister(17),
        trunc=notrunc(),
    )

    @test first_sample.configuration == repeated_sample.configuration
    @test first_sample.conditional_probabilities ≈ fill(0.5, number_sites)
    @test first_sample.log_probability ≈ -number_sites * log(2)
    @test first_sample.log_probability ≈
        sum(log, first_sample.conditional_probabilities)
    @test first_sample.compression_norm ≈ 1
    @test first_sample.maximum_probability_residual ≤ 1e-12
    @test first_sample.sector_list == collect(sectors(physical_space))
    @test norm(first_sample.state) ≈ 1

    for site in 1:number_sites
        tensor = first_sample.state[site]
        @test tensor isa MPSTensor{4}
        @test codomain(tensor)[2] == physical_space
        @test dim(domain(tensor)[1]) == 1
        if site < number_sites
            @test domain(tensor)[end] ==
                codomain(first_sample.state[site + 1])[1]
        end
    end
end

@testset "Z₂ complete-basis Born regression" begin
    physical_space = z2_born_space(0, 1)
    thermal_mpo = reference_thermal_mpo(physical_space)
    configurations = collect(Iterators.product(1:2, 1:2))
    direct_columns = [
        direct_column(thermal_mpo, physical_space, configuration)
        for configuration in configurations
    ]
    probabilities = getproperty.(direct_columns, :probability)

    @test norm(thermal_mpo) ≈ 1
    @test sum(probabilities) ≈ 1 atol = 1e-12
    @test length(unique(round.(probabilities; digits=10))) > 1

    purification = zeros(ComplexF64, length(configurations), length(configurations))
    normalized_columns = zeros(ComplexF64, size(purification))
    for (input_index, input_configuration) in enumerate(configurations)
        direct = direct_columns[input_index]
        for (output_index, output_configuration) in enumerate(configurations)
            output_parity = mod(sum(index - 1 for index in output_configuration), 2)
            input_parity = mod(sum(index - 1 for index in input_configuration), 2)
            if output_parity == input_parity
                basis = basis_state(
                    physical_space,
                    output_configuration,
                    input_configuration,
                )
                normalized_amplitude = inner(basis, direct.state)
                normalized_columns[output_index, input_index] = normalized_amplitude
                purification[output_index, input_index] =
                    sqrt(direct.probability) * normalized_amplitude
            end
        end
    end

    density_from_purification = purification * purification'
    density_from_born_columns = zeros(ComplexF64, size(density_from_purification))
    for input_index in eachindex(configurations)
        wavefunction = view(normalized_columns, :, input_index)
        density_from_born_columns .+= probabilities[input_index] .* (wavefunction * wavefunction')
    end
    density_trace = sum(
        density_from_purification[index, index]
        for index in axes(density_from_purification, 1)
    )
    @test density_trace ≈ 1 atol = 1e-12
    @test density_from_born_columns ≈ density_from_purification atol = 1e-12

    hamiltonian = dense_tfim_hamiltonian(configurations)
    direct_energy = real(sum(density_from_purification .* transpose(hamiltonian)))
    born_energy = sum(eachindex(configurations)) do input_index
        wavefunction = view(normalized_columns, :, input_index)
        sample_energy = real(sum(conj.(wavefunction) .* (hamiltonian * wavefunction)))
        probabilities[input_index] * sample_energy
    end
    @test born_energy ≈ direct_energy atol = 1e-12

    for (input_index, configuration) in enumerate(configurations)
        thresholds = configuration_thresholds(configuration, configurations, probabilities)
        sampled = z2_purification_born_sample(
            thermal_mpo,
            physical_space;
            rng=ThresholdRNG(thresholds),
            trunc=notrunc(),
        )
        direct = direct_columns[input_index]
        @test Tuple(sampled.configuration) == configuration
        @test exp(sampled.log_probability) ≈ direct.probability atol = 1e-12
        @test prod(sampled.conditional_probabilities) ≈ direct.probability atol = 1e-12
        @test abs(inner(sampled.state, direct.state)) ≈ 1 atol = 1e-12
    end
end
