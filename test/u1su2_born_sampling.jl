@testset "unconditioned U(1)×SU(2) Born backend" begin
    for thermal in (u1su2_reference_thermal_mpo(),)
        canonicalize!(thermal, 1)
        normalize!(thermal)
        distribution = zeros(5)
        for tuple in Iterators.product(1:3, 1:3)
            configuration = collect(tuple)
            born = u1su2_purification_born_project(thermal, configuration)
            direct = u1su2_direct_column(thermal, configuration)
            fixed = u1su2_fixed_n_purification_born_project(
                thermal, sum(configuration .- 1), configuration)
            @test born.born_probability ≈ direct.probability atol=1e-12
            @test born.born_probability ≈ fixed.born_probability * fixed.target_sector_weight atol=1e-12
            @test abs(inner(born.state, direct.state)) ≈ 1 atol=1e-12
            @test norm(born.state) ≈ 1 atol=1e-12
            for site in 1:2
                @test domain(born.state[site])[1] == StoMPS._born_sector_space(
                    u1su2_born_space(), born.sector_list[configuration[site]])
            end
            distribution[born.particle_number + 1] += born.born_probability
        end
        @test sum(distribution) ≈ 1 atol=1e-12
    end
    identity = identityMPO(2, u1su2_born_space())
    canonicalize!(identity, 1)
    normalize!(identity)
    distribution = zeros(5)
    for tuple in Iterators.product(1:3, 1:3)
        born = u1su2_purification_born_project(identity, collect(tuple))
        @test born.born_probability ≈ prod([1,2,1][i] for i in tuple) / 16 atol=1e-12
        distribution[born.particle_number + 1] += born.born_probability
    end
    @test distribution ≈ [1,4,6,4,1] / 16 atol=1e-12
    a = u1su2_purification_born_sample(identity; rng=MersenneTwister(12))
    b = u1su2_purification_born_sample(identity; rng=MersenneTwister(12))
    @test a.configuration == b.configuration
    @test_throws ArgumentError u1su2_purification_born_project(identity, [1,4])
end
