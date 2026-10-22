using StoMPS
using Test
using Random
using LinearAlgebra
using StoMPS.FiniteMPS
using StoMPS.FiniteMPS.TensorKit

isempty(ARGS) || throw(ArgumentError("usage: runtests.jl"))

@testset "StoMPS.jl" begin
    include("hubbard.jl")
    include("triangular_hubbard.jl")
end
