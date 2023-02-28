using LinearAlgebra: dot, I
using Cones
using Cones: SOCone, GeneralSOCone
using Test

@testset "SOCone membership" begin
    K1 = SOCone(3, 1.0)
    @test [1.0, 0.0, 0.0] ∈ K1
    @test [0.0, 1.0, 0.0] ∉ K1
end

@testset "SOCone on_boundary" begin
    K1 = SOCone(3, 1.0)
    @test !(Cones.on_boundary([1.0, 0.0, 0.0], K1))
    @test Cones.on_boundary([0.0, 0.0, 0.0], K1)
    @test Cones.on_boundary([1.0, 1/sqrt(2), 1/sqrt(2)], K1)
end

@testset "SOCone move_to_boundary" begin
    K1 = SOCone(3, 2.0)
    p1 = Cones.move_to_boundary([1,0, 0.0, 0.0], K1) # Will be pushed to origin
    p2 = Cones.move_to_boundary([0.0, 1.0, 1.0], K1) # Point will be pushed upwards from outside of cone to boundary
    p3 = Cones.move_to_boundary([4.0, 0.5, 1.0], K1) # Point inside cone will be pushed downwards to boundary
    @test Cones.on_boundary(p1, K1)
    @test Cones.on_boundary(p2, K1)
    @test Cones.on_boundary(p3, K1)
end

@testset "SOCone - duality" begin
    K1 = SOCone(4, 0.5)
    K2 = Cones.dual(K1)
    for i in 1:10
        p = Cones.move_to_boundary(rand(K1.d), K1) # Create point on cone boundary
        r = copy(p)
        r[2:K1.d] .*= -1 # reflection with respect to cone axis
        q = Cones.move_to_boundary(r, K2) # Point now on boundary of dual cone on opposity side
        @test dot(p, q) ≈ 0.0 atol = sqrt(eps()) # Constructed points should be orthogonal

        # Point "nudged downwards" should not be outside of dual cone
        # and have a negative dot product with p1
        # Note that this only works if the initial cone contained in dual
        # i.e. α <= 1
        t = q .- 1e-3 
        @test dot(p, t) < 0.0
    end
end

@testset "GeneralSOCone membership - standard SOCone" begin
    K1 = GeneralSOCone([1, 0, 0], 1.0)
    @test [1.0, 0.0, 0.0] ∈ K1
    @test [0.0, 1.0, 0.0] ∉ K1

    K2 = GeneralSOCone([1.0/sqrt(2), -1.0/sqrt(2)], 1.0) # SOCone rotated 45 degrees clockwise
    @test [0.0, 1.0] ∉ K2 # not in original and not transformed
    @test [1.0, 0.0] ∈ K2 # in original and transformed
    @test [0.0, -1.0] ∈ K2 # not in original but in transformed
    @test [1.0, 1.0] ∉ K2 # in original but not in transformed
end


@testset "TransformedCone membership" begin
    C = SOCone(2, 1.0)
    U = (I/sqrt(2)) * [[1.0, -1] [1, 1]]
    K = TransformedCone(C, U) # SOCone rotated 90 degrees clockwise
    @test [0.0, 1.0] ∉ K # not in original and not transformed
    @test [1.0, 0.0] ∈ K # in original and transformed
    @test [0.0, -1.0] ∈ K # not in original but in transformed
    @test [1.0, 1.0] ∉ K # in original but not in transformed
end