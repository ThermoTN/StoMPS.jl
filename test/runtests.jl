using StoMPS
using Test

@testset "StoMPS.jl" begin
    include("estimators.jl")
    include("z2_born_sampling.jl")
    include("z2_born_ed.jl")
    include("su2su2_born_sampling.jl")
    include("u1su2_fixed_n_born_sampling.jl")
    include("u1su2_born_sampling.jl")
    include("su2su2_resampling.jl")
    include("resampling_mixture.jl")
    include("adaptive_resampling.jl")
    include("su2su2_born_ed.jl")
    include("born_continuation.jl")
    include("hubbard.jl")
    include("triangular_hubbard.jl")
    include("workflow.jl")
    include("dynamics.jl")
end
