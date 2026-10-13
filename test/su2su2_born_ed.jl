function dense_creation_operator(mode, number_modes)
    dimension = 1 << number_modes
    result = zeros(Float64, dimension, dimension)
    occupied_bit = 1 << (mode - 1)
    lower_modes = occupied_bit - 1
    for state in 0:(dimension - 1)
        iszero(state & occupied_bit) || continue
        sign = isodd(count_ones(state & lower_modes)) ? -1.0 : 1.0
        result[(state | occupied_bit) + 1, state + 1] = sign
    end
    return result
end

function dense_two_site_hubbard(; hopping=1.0, interaction=4.0)
    number_sites = 2
    number_modes = 2 * number_sites
    creators = [dense_creation_operator(mode, number_modes) for mode in 1:number_modes]
    hamiltonian = zeros(Float64, 1 << number_modes, 1 << number_modes)
    for spin in 0:1
        first = creators[spin + 1]
        second = creators[2 + spin + 1]
        hamiltonian .-= hopping .* (first' * second + second' * first)
    end
    chemical_potential = interaction / 2
    for site in 1:number_sites
        up = creators[2 * site - 1]' * creators[2 * site - 1]
        down = creators[2 * site]' * creators[2 * site]
        hamiltonian .+= interaction .* (up * down)
        hamiltonian .-= chemical_potential .* (up + down)
    end
    return hamiltonian
end

function dense_su2su2_sector_snapshot(configuration, sector_list)
    spin_sector = first(sector_list)
    local_states = map(configuration) do sector_index
        sector_list[sector_index] == spin_sector ? (1, 2) : (0, 3)
    end
    basis_indices = Int[]
    for local_configuration in Iterators.product(local_states...)
        fock_state = sum(
            local_state << (2 * (site - 1))
            for (site, local_state) in enumerate(local_configuration)
        )
        push!(basis_indices, fock_state + 1)
    end
    snapshot = zeros(Float64, 16, length(basis_indices))
    for (column, basis_index) in enumerate(basis_indices)
        snapshot[basis_index, column] = 1 / sqrt(length(basis_indices))
    end
    return snapshot
end

@testset "SO(4) Born energy per site agrees with Hubbard ED" begin
    number_sites = 2
    target_beta = 0.7
    physical_space = su2su2_born_space()
    sector_list = collect(sectors(physical_space))
    configurations = collect(Iterators.product(1:2, 1:2))
    hamiltonian = dense_two_site_hubbard()
    half_evolution = exp(-target_beta * hamiltonian / 2)

    identity_purification = identityMPO(number_sites, physical_space)
    normalize!(identity_purification)
    exact_born_probability = 1 / length(configurations)
    sample_energies = zeros(Float64, length(configurations))
    deterministic_log_weights = similar(sample_energies)

    for (sample_index, configuration) in enumerate(configurations)
        thresholds = configuration_thresholds(
            configuration,
            configurations,
            fill(exact_born_probability, length(configurations)),
        )
        born = su2su2_purification_born_sample(
            identity_purification,
            physical_space;
            rng=ThresholdRNG(thresholds),
            trunc=notrunc(),
        )
        @test Tuple(born.configuration) == configuration
        @test exp(born.log_probability) ≈ exact_born_probability atol = 1e-12

        initial_snapshot = dense_su2su2_sector_snapshot(configuration, sector_list)
        @test norm(initial_snapshot) ≈ 1 atol = 1e-12
        evolved_snapshot = half_evolution * initial_snapshot
        evolved_norm = norm(evolved_snapshot)
        deterministic_log_weights[sample_index] =
            log(exact_born_probability) + 2 * log(evolved_norm)
        sample_energies[sample_index] = real(
            tr(evolved_snapshot' * hamiltonian * evolved_snapshot),
        ) / abs2(evolved_norm)
    end

    born_energy_per_site = logweighted_mean(
        sample_energies,
        deterministic_log_weights,
    ) / number_sites
    gibbs_operator = exp(-target_beta * hamiltonian)
    ed_energy_per_site = real(tr(hamiltonian * gibbs_operator)) /
        real(tr(gibbs_operator)) / number_sites
    absolute_error_per_site = abs(born_energy_per_site - ed_energy_per_site)
    @test absolute_error_per_site <= 1e-9
end
