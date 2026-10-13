using LinearAlgebra

"""Dense two-site oracle for Appendix-B posterior mixture sampling.

The two local values are the complete eta/spin sector labels; no magnetic
quantum number is represented in this oracle.
"""
function dense_configuration_index(configuration)
    return 2 * (configuration[1] - 1) + configuration[2]
end

function dense_member_prefix_probability(state, prefix)
    probability = 0.0
    for configuration in Iterators.product(1:2, 1:2)
        all(configuration[site] == prefix[site] for site in eachindex(prefix)) ||
            continue
        probability += abs2(state[dense_configuration_index(configuration)])
    end
    return probability
end

function dense_empirical_probability(states, log_weights, configuration)
    weights = exp.(log_weights .- maximum(log_weights))
    return sum(weights[index] * abs2(states[index][dense_configuration_index(configuration)])
        for index in eachindex(states)) / sum(weights)
end

function dense_empirical_prefix_probability(states, log_weights, prefix)
    weights = exp.(log_weights .- maximum(log_weights))
    return sum(weights[member] * dense_member_prefix_probability(states[member], prefix)
        for member in eachindex(states)) / sum(weights)
end

function posterior_sweep_probability(states, log_weights, configuration)
    posterior = copy(log_weights)
    conditional_probabilities = Float64[]
    prefix = Int[]
    for site in eachindex(configuration)
        member_likelihoods = zeros(Float64, length(states), 2)
        for member in eachindex(states)
            isfinite(posterior[member]) || continue
            prefix_mass = dense_member_prefix_probability(states[member], prefix)
            prefix_mass == 0 && continue
            for sector in 1:2
                member_likelihoods[member, sector] =
                    dense_member_prefix_probability(states[member], [prefix; sector]) /
                    prefix_mass
            end
        end
        shift = maximum(posterior)
        normalized = exp.(posterior .- shift)
        normalized ./= sum(normalized)
        mixture = vec(transpose(normalized) * member_likelihoods)
        push!(conditional_probabilities, mixture[configuration[site]])
        likelihoods = member_likelihoods[:, configuration[site]]
        for member in eachindex(posterior)
            posterior[member] = likelihoods[member] == 0 ? -Inf :
                posterior[member] + log(likelihoods[member])
        end
        push!(prefix, configuration[site])
    end
    return (; probability=prod(conditional_probabilities), conditional_probabilities)
end

function fixed_weight_sweep_probability(states, log_weights, configuration)
    weights = exp.(log_weights .- maximum(log_weights))
    weights ./= sum(weights)
    prefix = Int[]
    probability = 1.0
    for site in eachindex(configuration)
        mixture_probability = sum(weights[member] * begin
            prefix_mass = dense_member_prefix_probability(states[member], prefix)
            prefix_mass == 0 ? 0.0 :
                dense_member_prefix_probability(states[member], [prefix; configuration[site]]) /
                prefix_mass
        end for member in eachindex(states))
        probability *= mixture_probability
        push!(prefix, configuration[site])
    end
    return probability
end

@testset "Explicit posterior empirical-mixture oracle" begin
    states = [
        normalize(ComplexF64[0.8, 0.1, 0.1, 0.4]),
        normalize(ComplexF64[0.1, 0.7, 0.6, 0.2]),
    ]
    log_weights = [log(1.0), log(3.0)]
    configurations = collect(Iterators.product(1:2, 1:2))
    exact = [dense_empirical_probability(states, log_weights, configuration)
        for configuration in configurations]
    sweeps = [posterior_sweep_probability(states, log_weights, configuration)
        for configuration in configurations]
    swept = getfield.(sweeps, :probability)
    @test swept ≈ exact atol = 1e-12
    @test sum(swept) ≈ 1 atol = 1e-12
    # Each local conditional is checked against the corresponding prefix of
    # the dense empirical mixed state, not merely against the final outcome.
    for (configuration, sweep) in zip(configurations, sweeps)
        running_probability = 1.0
        for prefix_length in 1:2
            running_probability *= sweep.conditional_probabilities[prefix_length]
            @test running_probability ≈ dense_empirical_prefix_probability(
                states,
                log_weights,
                collect(configuration[1:prefix_length]),
            ) atol = 1e-12
        end
    end
    # This proves that keeping the initial member weights after a prefix is not
    # an equivalent implementation of the Appendix-B conditional mixture.
    wrong = [fixed_weight_sweep_probability(states, log_weights, configuration)
        for configuration in configurations]
    @test maximum(abs.(wrong .- exact)) > 1e-3

    zero_branch_states = [
        ComplexF64[1, 0, 0, 0],
        normalize(ComplexF64[0, 1, 1, 1]),
    ]
    zero_sweep = posterior_sweep_probability(
        zero_branch_states,
        [0.0, 0.0],
        (2, 1),
    )
    @test zero_sweep.probability ≈ dense_empirical_probability(
        zero_branch_states,
        [0.0, 0.0],
        (2, 1),
    ) atol = 1e-12
end

@testset "Regenerated proposal correction and absolute logZ" begin
    p_beta0 = [0.2, 0.3, 0.5]
    proposal = [0.5, 0.25, 0.25]
    full_norm = [2.0, 3.0, 4.0]
    energies = [-3.0, -2.0, -1.0]
    log_importance = log.(p_beta0) .- log.(proposal) .+ log.(full_norm)

    @test sum(proposal .* exp.(log_importance)) ≈ sum(p_beta0 .* full_norm)
    exact_energy = sum(p_beta0 .* full_norm .* energies) / sum(p_beta0 .* full_norm)
    draws = vcat(fill(1, 50), fill(2, 25), fill(3, 25))
    sampled_log_weights = log_importance[draws]
    @test logweighted_mean(energies[draws], sampled_log_weights) ≈ exact_energy
    logZ_beta0 = 2.75
    @test estimate_logpartition(logZ_beta0, sampled_log_weights) ≈
        logZ_beta0 + log(sum(p_beta0 .* full_norm))

    # Omitting p_beta0 / proposal is a different, biased estimator.
    @test logweighted_mean(energies[draws], log.(full_norm[draws])) != exact_energy

    target_proposal = p_beta0 .* full_norm ./ sum(p_beta0 .* full_norm)
    ideal_log_importance = log.(p_beta0) .- log.(target_proposal) .+ log.(full_norm)
    @test sampling_efficiency(ideal_log_importance) ≈ 1.0 atol = 1e-14
end
