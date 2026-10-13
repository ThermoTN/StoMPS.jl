function u1su2_reference_thermal_mpo()
    physical_space = u1su2_born_space()
    thermal = identityMPO(2, physical_space)
    sector_list = StoMPS._u1su2_sector_list()
    weights = ((0.8, 1.1, 1.4), (1.3, 0.9, 0.7))
    for site in eachindex(thermal), sector_index in eachindex(sector_list)
        block(thermal[site].A, sector_list[sector_index]) .*= weights[site][sector_index]
    end
    canonicalize!(thermal, length(thermal); trunc=notrunc())
    canonicalize!(thermal, 1; trunc=notrunc())
    normalize!(thermal)
    return thermal
end

function u1su2_direct_column(thermal, configuration)
    physical_space = u1su2_born_space()
    sector_list = StoMPS._u1su2_sector_list()
    tensors = AbstractMPSTensor[]
    for site in eachindex(configuration)
        push!(tensors, StoMPS._u1su2_projected_tensor(
            thermal[site], physical_space, sector_list[configuration[site]],
        ))
    end
    state = MPS(tensors, [1, length(tensors)], FiniteMPS.coef(thermal))
    canonicalize!(state, length(state); trunc=notrunc())
    canonicalize!(state, 1; trunc=notrunc())
    probability = abs2(norm(state))
    probability > 0 && normalize!(state)
    return (; state, probability)
end

function u1su2_full_chain_sector_probability(state, site, sector_index)
    physical_space = u1su2_born_space()
    sector = StoMPS._u1su2_sector_list()[sector_index]
    marker_space = StoMPS._born_sector_space(physical_space, sector)
    embedding = isometry(physical_space, marker_space)
    projector = embedding * embedding'
    projected = deepcopy(state)
    center = projected[site].A
    @tensor tensor[left physical; marker right] :=
        projector[physical; input] * center[left input; marker right]
    projected[site] = MPSTensor(tensor)
    return real(inner(projected, projected))
end

@testset "fixed-N U(1)×SU(2) Born backend" begin
    physical_space = u1su2_born_space()
    sector_list = StoMPS._u1su2_sector_list()
    @test [Int(sector[1].charge) for sector in sector_list] == [-1, 0, 1]
    @test [dim(StoMPS._born_sector_space(physical_space, sector))
           for sector in sector_list] == [1, 2, 1]
    thermal = u1su2_reference_thermal_mpo()

    all_configurations = collect(Iterators.product(1:3, 1:3))
    direct = Dict(configuration => u1su2_direct_column(thermal, configuration)
                  for configuration in all_configurations)
    for target_N in 0:4
        configurations = filter(configuration ->
            sum(index -> index - 1, configuration) == target_N,
            all_configurations,
        )
        sector_weight = sum(direct[configuration].probability
                            for configuration in configurations)
        conditional_sum = 0.0
        for configuration in configurations
            projected = u1su2_fixed_n_purification_born_project(
                thermal, target_N, collect(configuration); trunc=notrunc(),
            )
            expected = direct[configuration].probability / sector_weight
            conditional_sum += projected.born_probability
            @test projected.occupancies == collect(configuration) .- 1
            @test projected.marker_dimensions == [1, 2, 1]
            @test projected.target_sector_weight ≈ sector_weight atol=1e-12
            @test projected.born_probability ≈ expected atol=1e-12
            @test exp(projected.log_probability) ≈ expected atol=1e-12
            @test abs(inner(projected.state, direct[configuration].state)) ≈ 1 atol=1e-12
            @test norm(projected.state) ≈ 1 atol=1e-12
        end
        @test conditional_sum ≈ 1 atol=1e-12
    end

    sample = u1su2_fixed_n_purification_born_sample(
        thermal, 2; rng=MersenneTwister(20260902), trunc=notrunc(),
    )
    repeated = u1su2_fixed_n_purification_born_sample(
        thermal, 2; rng=MersenneTwister(20260902), trunc=notrunc(),
    )
    @test sample.configuration == repeated.configuration
    @test sum(sample.occupancies) == 2

    branch = u1su2_resampling_branch(sample.state)
    for site in 1:length(sample.state)
        rdm_local = u1su2_resampling_local_rdm(branch, site)
        @test sum(rdm_local.sector_probabilities) ≈ 1 atol=1e-12
        @test rdm_local.marker_dimensions == [1, 2, 1]
        @test rdm_local.sector_probabilities ≈ [
            u1su2_full_chain_sector_probability(branch.state, site, sector_index)
            for sector_index in 1:3
        ] atol=1e-12
        selected = findfirst(>(1e-14), rdm_local.sector_probabilities)
        result = u1su2_resampling_project!(branch, site, selected)
        @test result.active
    end

    @test_throws ArgumentError u1su2_fixed_n_purification_born_project(
        thermal, 2, [1, 1],
    )
    @test_throws ArgumentError u1su2_fixed_n_purification_born_sample(
        thermal, 5; rng=MersenneTwister(1),
    )
end
