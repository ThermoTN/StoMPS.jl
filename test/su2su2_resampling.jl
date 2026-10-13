function full_chain_sector_probability(state, physical_space, site, sector_index)
    sector = collect(sectors(physical_space))[sector_index]
    marker_space = StoMPS._born_sector_space(physical_space, sector)
    embedding = isometry(physical_space, marker_space)
    projector = embedding * embedding'
    projected_state = deepcopy(state)
    center = projected_state[site].A
    @tensor projected[left physical; marker right] :=
        projector[physical; input] * center[left input; marker right]
    projected_state[site] = MPSTensor(projected)
    # `norm(MPS)` may use the stored center/isometry metadata; `inner` forces
    # the independent whole-chain partial trace used as this test oracle.
    return real(inner(projected_state, projected_state))
end

@testset "SO(4) annealed-resampling tensor kernel" begin
    physical_space = su2su2_born_space()
    thermal_mpo = su2su2_reference_thermal_mpo(physical_space)
    configuration = [1, 2]
    fixed = su2su2_purification_born_project(
        thermal_mpo,
        physical_space,
        configuration;
        trunc=notrunc(),
    )
    direct = su2su2_direct_sector_column(thermal_mpo, physical_space, Tuple(configuration))
    @test fixed.configuration == configuration
    @test exp(fixed.log_probability) ≈ direct.probability atol = 1e-12
    @test abs(inner(fixed.state, direct.state)) ≈ 1 atol = 1e-12

    branch = su2su2_resampling_branch(fixed.state, physical_space)
    local_one = su2su2_resampling_local_rdm(branch, 1)
    @test local_one.trace_residual ≤ 1e-12
    @test local_one.hermiticity_residual ≤ 1e-12
    @test local_one.probability_residual ≤ 1e-12
    @test sum(local_one.sector_probabilities) ≈ 1 atol = 1e-12
    @test local_one.sector_probabilities ≈ [
        full_chain_sector_probability(fixed.state, physical_space, 1, index)
        for index in 1:2
    ] atol = 1e-12

    selected = findfirst(>(0), local_one.sector_probabilities)
    projection = su2su2_resampling_project!(branch, 1, selected)
    @test projection.likelihood ≈ local_one.sector_probabilities[selected] atol = 1e-12
    @test projection.active
    @test Center(branch.state) == [2, 2]
    local_two = su2su2_resampling_local_rdm(branch, 2)
    selected_two = findfirst(>(0), local_two.sector_probabilities)
    su2su2_resampling_project!(branch, 2, selected_two)
    @test norm(branch.state) ≈ 1 atol=1e-12

    off_center = deepcopy(fixed.state)
    canonicalize!(off_center, 2; trunc=notrunc())
    @test Center(off_center) == [2, 2]
    fallback = su2su2_resampling_branch(off_center, physical_space)
    @test Center(fallback.state) == [1, 1]
    @test su2su2_resampling_local_rdm(fallback, 1).sector_probabilities ≈
        local_one.sector_probabilities atol = 1e-12
end
