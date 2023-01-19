using Test
using LinearAlgebra: I
using Cones

@testset "TransformedCone membership" begin
    C = SOCone(2, 1.0)
    U = (I/sqrt(2)) * [[1.0, -1] [1, 1]]
    K = TransformedCone(C, U) # SOCone rotated 45 degrees clockwise
    @test [0.0, 1.0] ∉ K # not in original and not transformed
    @test [1.0, 0.0] ∈ K # in original and transformed
    @test [0.0, -1.0] ∈ K # not in original but in transformed
    @test [1.0, 1.0] ∉ K # in original but not in transformed
end