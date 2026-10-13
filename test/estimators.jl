@testset "Finite-temperature estimators" begin
    values = [1.0, 3.0]
    log_weights = [1000.0, 1000.0 + log(3)]
    @test logweighted_mean(values, log_weights) ≈ 2.5
    @test logweighted_mean(values, log_weights .- 1234) ≈ 2.5
    @test estimate_logpartition(2.0, [0.0, log(3)]) ≈ 2.0 + log(2)

    @test sampling_efficiency(zeros(4)) == 1.0
    expected_efficiency = (1 + 2 + 4)^2 / (3 * (1 + 4 + 16))
    @test sampling_efficiency([0.0, log(2), log(4)]) ≈ expected_efficiency
    @test sampling_efficiency([0.0, log(2), log(4)] .+ 1000) ≈
        expected_efficiency
    @test sampling_efficiency([0.0, -1_000.0, -1_000.0]) ≈ 1 / 3

    @test_throws ArgumentError logweighted_mean([Inf], [0.0])

    # Independent two-point delta-method oracle and virtual-copy equivalence.
    errors=StoMPS.conditional_weight_errors(log.([1.,3.]),[2,1])
    expanded=StoMPS.conditional_weight_errors(log.([1.,1.,3.]),ones(Int,3))
    @test errors.lnZ_error ≈ 0.4
    @test errors.eta_error ≈ 20/363
    @test errors == expanded
    @test StoMPS.conditional_weight_errors(zeros(2),[3,2]) == (lnZ_error=0.0,eta_error=0.0)
end
