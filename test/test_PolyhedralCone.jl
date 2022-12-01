# Test whether PolyhedralCone projects points onto same locations as equivalent FiniteCones

using Distributions
using Random

Random.seed!(1)

print("\tTest 1 - dim == 2...")

@testset "projection onto polyhedral cone - equivalence with finite cone" begin
    Z = Normal()
    μ = rand(Z, 2)
    A = rand(Z, 2, 2)
    Σ = A'A
    dist = MvNormal(μ, Σ)

    K_f = FiniteCone([1 0; 0 1])
    K_p = PolyhedralCone([1 0; 0 1])

    for i in 1:1000
        x = rand(dist)
        @test project(K_f, x) ≈ project(K_p, x) atol=1e-4
    end
end

@testset "membership of polyhedral cone - simple examples" begin
    K1 = PolyhedralCone([2 -1; -1 2]) # (y <= 2x, y >= (1/2)x )
    @test [1.0,1] ∈ K1
    @test [1.0,0] ∉ K1
    @test [0,1] ∉ K1
    @test [1,2] ∈ K1
    @test [2,1] ∈ K1
    @test [-1, -1] ∉ K1
end

