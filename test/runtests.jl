using GeneralizedMultiparticleMie
using GSL
using Random
using Test
using WignerD

Random.seed!(42)

@testset "GeneralizedMultiparticleMie.jl" begin
    include("test_special_functions.jl")
end
