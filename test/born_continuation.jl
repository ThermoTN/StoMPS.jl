using Random
using StoMPS.FiniteMPS

@testset "Born sample continuation" begin
    physical_space = z2_born_space(0, 1)
    thermal_mpo = identityMPO(2, physical_space)
    normalize!(thermal_mpo)
    born = z2_purification_born_sample(
        thermal_mpo,
        physical_space;
        rng=MersenneTwister(11),
    )

    beta_grid = [0.5, 0.75, 1.0]
    scales = [1.0, 2.0, 0.5]
    time_steps = Float64[]
    random_draws = Float64[]
    trajectory = continue_born_sample!(
        born.state,
        beta_grid;
        rng=MersenneTwister(29),
        step! = function (sample, time_step, beta_index, rng)
            push!(time_steps, time_step)
            push!(random_draws, rand(rng))
            rmul!(sample, scales[beta_index])
        end,
        observe=(sample, beta, beta_index) -> beta + beta_index - beta_index,
    )

    reference_rng = MersenneTwister(29)
    @test trajectory.beta == beta_grid
    @test trajectory.observable == beta_grid
    @test trajectory.relative_log_weight ≈ [0.0, 2 * log(2), 0.0]
    @test time_steps ≈ [-0.125, -0.125]
    @test random_draws == [rand(reference_rng), rand(reference_rng)]
    @test norm(born.state) ≈ 1

    no_op_step!(sample, time_step, beta_index, rng) = nothing
    observe_norm(sample, beta, beta_index) = norm(sample)
    unnormalized = deepcopy(born.state)
    rmul!(unnormalized, 2)
    @test_throws ArgumentError continue_born_sample!(
        unnormalized,
        [0.5];
        rng=MersenneTwister(1),
        step! = no_op_step!,
        observe=observe_norm,
    )
    @test_throws ArgumentError continue_born_sample!(
        deepcopy(born.state),
        [0.5];
        rng=MersenneTwister(1),
        step! = no_op_step!,
        observe=(sample, beta, beta_index) -> Inf,
    )
end
