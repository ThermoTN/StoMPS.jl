using LinearAlgebra

function dense_snapshot(
    born_state,
    physical_space,
    configurations,
    input_configuration,
)
    wavefunction = zeros(ComplexF64, length(configurations))
    input_parity = mod(sum(index - 1 for index in input_configuration), 2)
    for (output_index, output_configuration) in enumerate(configurations)
        output_parity = mod(sum(index - 1 for index in output_configuration), 2)
        if output_parity == input_parity
            basis = basis_state(
                physical_space,
                output_configuration,
                input_configuration,
            )
            wavefunction[output_index] = inner(basis, born_state)
        end
    end
    return wavefunction
end

@testset "Z₂ Born per-site energy agrees with ED" begin
    number_sites = 2
    target_beta = 0.7
    configurations = collect(Iterators.product(1:2, 1:2))
    hamiltonian = dense_tfim_hamiltonian(
        configurations;
        coupling=0.8,
        field=0.35,
    )

    physical_space = z2_born_space(0, 1)
    identity_purification = identityMPO(number_sites, physical_space)
    logZ0 = 2 * log(norm(identity_purification))
    normalize!(identity_purification)

    exact_born_probabilities = fill(1 / length(configurations), length(configurations))
    sample_energies = zeros(Float64, length(configurations))
    relative_log_weights = similar(sample_energies)
    half_evolution = exp(-target_beta * hamiltonian / 2)

    for (sample_index, configuration) in enumerate(configurations)
        thresholds = configuration_thresholds(
            configuration,
            configurations,
            exact_born_probabilities,
        )
        born = z2_purification_born_sample(
            identity_purification,
            physical_space;
            rng=ThresholdRNG(thresholds),
            trunc=notrunc(),
        )
        @test Tuple(born.configuration) == configuration

        initial_state = dense_snapshot(
            born.state,
            physical_space,
            configurations,
            configuration,
        )
        @test norm(initial_state) ≈ 1 atol = 1e-12

        evolved_state = half_evolution * initial_state
        evolved_norm = norm(evolved_state)
        relative_log_weights[sample_index] = 2 * log(evolved_norm)
        sample_energies[sample_index] = real(
            dot(evolved_state, hamiltonian * evolved_state),
        ) / abs2(evolved_norm)
    end

    born_energy_per_site = logweighted_mean(
        sample_energies,
        relative_log_weights,
    ) / number_sites

    gibbs_operator = exp(-target_beta * hamiltonian)
    partition_function = real(tr(gibbs_operator))
    ed_energy_per_site = real(tr(hamiltonian * gibbs_operator)) /
        partition_function / number_sites
    absolute_error_per_site = abs(born_energy_per_site - ed_energy_per_site)

    # Complete enumeration has no Monte Carlo uncertainty.  The disjunction is
    # kept in the production acceptance form for future sampled regressions.
    energy_errorbar_per_site = 0.0
    @test absolute_error_per_site <= 1e-9 ||
        absolute_error_per_site <= energy_errorbar_per_site

    estimated_logZ = estimate_logpartition(logZ0, relative_log_weights)
    @test estimated_logZ ≈ log(partition_function) atol = 1e-12
end
