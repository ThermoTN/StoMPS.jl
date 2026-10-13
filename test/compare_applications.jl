using Test, StoMPS
include("workflow_helpers.jl")
reference=load_result(joinpath(ARGS[1],"results","result.jld2"))
@testset "MPI rank-independent trajectories" begin
    for directory in ARGS[2:end]
        other=load_result(joinpath(directory,"results","result.jld2"))
        compare_curves(reference["observables"]["curve"],other["observables"]["curve"])
    end
end
