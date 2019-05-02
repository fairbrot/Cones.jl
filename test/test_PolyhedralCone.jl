# Test whether PolyhedralCone projects points onto same locations as equivalent FiniteCones

using Distributions
using Random

Random.seed!(1)

print("\tTest 1 - dim == 2...")
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

println("success!")
