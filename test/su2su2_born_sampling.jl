function su2su2_marker_space(physical_space, sector)
    return StoMPS._born_sector_space(physical_space, sector)
end

function su2su2_reference_thermal_mpo(physical_space)
    spin_sector, charge_sector = collect(sectors(physical_space))
    thermal_mpo = identityMPO(2, physical_space)
    local_weights = ((0.7, 1.3), (1.1, 0.6))
    for (site, (spin_weight, charge_weight)) in enumerate(local_weights)
        block(thermal_mpo[site].A, spin_sector) .*= spin_weight
        block(thermal_mpo[site].A, charge_sector) .*= charge_weight
    end
    canonicalize!(thermal_mpo, length(thermal_mpo); trunc=notrunc())
    canonicalize!(thermal_mpo, 1; trunc=notrunc())
    normalize!(thermal_mpo)
    return thermal_mpo
end

function su2su2_direct_sector_column(thermal_mpo, physical_space, configuration)
    sector_list = collect(sectors(physical_space))
    tensors = AbstractMPSTensor[]
    for site in eachindex(configuration)
        sector = sector_list[configuration[site]]
        marker_space = su2su2_marker_space(physical_space, sector)
        embedding = isometry(physical_space, marker_space)
        local_tensor = thermal_mpo[site].A
        @tensor projected[left physical; marker right] :=
            local_tensor[left physical; input right] * embedding[input; marker]
        push!(tensors, MPSTensor(projected))
    end
    column = MPS(tensors, [1, length(tensors)], FiniteMPS.coef(thermal_mpo))
    canonicalize!(column, length(column); trunc=notrunc())
    canonicalize!(column, 1; trunc=notrunc())
    probability = abs2(norm(column))
    probability > 0 && normalize!(column)
    return (; state=column, probability)
end

@testset "SO(4) degenerate purification Born sampling" begin
    physical_space = su2su2_born_space()
    sector_list = collect(sectors(physical_space))
    @test string.(sector_list) == [
        "Irrep[SU₂ × SU₂](0, 1/2)",
        "Irrep[SU₂ × SU₂](1/2, 0)",
    ]
    @test all(sector -> dim(physical_space, sector) == 1, sector_list)
    @test all(sector -> dim(su2su2_marker_space(physical_space, sector)) == 2, sector_list)

    projectors = map(sector_list) do sector
        embedding = isometry(physical_space, su2su2_marker_space(physical_space, sector))
        embedding * embedding'
    end
    @test projectors[1] + projectors[2] ≈ id(physical_space)
    @test norm(projectors[1] * projectors[2]) ≤ 1e-12
    @test StoMPS._draw_born_probability_index(
        [0.0, 1.0],
        ThresholdRNG([0.0]),
    ) == 2

    unnormalized = identityMPO(2, physical_space)
    @test_throws ArgumentError su2su2_purification_born_sample(
        unnormalized,
        physical_space;
        rng=MersenneTwister(1),
    )
    wrong_space = identityMPO(2, SU2Spin.pspace)
    normalize!(wrong_space)
    @test_throws ArgumentError su2su2_purification_born_sample(
        wrong_space,
        SU2Spin.pspace;
        rng=MersenneTwister(1),
    )

    thermal_mpo = su2su2_reference_thermal_mpo(physical_space)
    configurations = collect(Iterators.product(1:2, 1:2))
    direct_columns = [
        su2su2_direct_sector_column(thermal_mpo, physical_space, configuration)
        for configuration in configurations
    ]
    probabilities = getproperty.(direct_columns, :probability)
    @test sum(probabilities) ≈ 1 atol = 1e-12
    @test length(unique(round.(probabilities; digits=10))) > 1

    for (configuration_index, configuration) in enumerate(configurations)
        thresholds = configuration_thresholds(configuration, configurations, probabilities)
        sampled = su2su2_purification_born_sample(
            thermal_mpo,
            physical_space;
            rng=ThresholdRNG(thresholds),
            trunc=notrunc(),
        )
        direct = direct_columns[configuration_index]
        @test Tuple(sampled.configuration) == configuration
        @test sampled.marker_dimensions == [2, 2]
        @test exp(sampled.log_probability) ≈ direct.probability atol = 1e-12
        @test prod(sampled.conditional_probabilities) ≈ direct.probability atol = 1e-12
        @test abs(inner(sampled.state, direct.state)) ≈ 1 atol = 1e-12
        @test sampled.maximum_probability_residual ≤ 1e-12
        @test norm(sampled.state) ≈ 1 atol = 1e-12
        for site in eachindex(configuration)
            @test domain(sampled.state[site])[1] ==
                su2su2_marker_space(physical_space, sector_list[configuration[site]])
        end
    end
end
