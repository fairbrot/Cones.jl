# Test whether PolyhedralCone projects points onto same locations as equivalent FiniteCones

using Distributions

srand(1)

print("\tTest 1 - dim == 2...")
Z = Normal()
μ = rand(Z, 2)
A = rand(Z, 2, 2)
Σ = A'A
dist = MvNormal(μ, Σ)

K_f = FiniteCone(eye(2))
K_p = PolyhedralCone(eye(2))

for i in 1:1000
    x = rand(dist)
    @test project(K_f, x) ≈ project(K_p, x) atol=1e-4
end

println("success!")
