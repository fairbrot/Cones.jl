using Random
using LinearAlgebra: rank

function test_gram_schmidt(A::Matrix{T}) where T<:Integer
    B = Cones.gram_schmidt(A)
    @test Cones.check_column_span_same(A, B)
    @test Cones.check_columns_orthogonal(B)
end

Random.seed!(1)

println("\t...gram_schmidt")

N = 10         # Trials
m, n = 5, 3
for counter in 0:(N-1)
    # Randomly generate full column rank matrix
    A = rand(Int, m, n) .% 5
    if rank(A) < n continue end
    test_gram_schmidt(A)
end

# Test normalise_rays
# 
function test_normalise_rays(rays::Matrix{T}, basis::Matrix{T}, norm_rays::Matrix{T}) where T<:Integer
    orth_basis = Cones.gram_schmidt(basis)
    norm_rays_out = Cones.normalise_rays(rays, basis)
    @test Cones.check_columns_same(norm_rays, norm_rays_out)
end

println("\t...normalise_rays")

# 0D lineality space
rays = [[2, 3, 1] [-1, 3, 2]]
basis = Array{Int}(undef, 3, 0)
norm_rays = rays
test_normalise_rays(rays, basis, norm_rays)

# 1D lineality space
rays = [[2, 0, 2] [0, 1, 1]]
basis= reshape([0, 0, 1], 3, 1)
norm_rays = [[1, 0, 0] [0, 1, 0]]
test_normalise_rays(rays, basis, norm_rays)

rays = [[0, 2, 1] [3, 1, 1]]
basis = reshape([-3, -3, 0], 3, 1)
norm_rays = [[-1, 1, 1] [1, -1, 1]]
test_normalise_rays(rays, basis, norm_rays)
         
# 2D lineality space
rays = reshape([1, 1, 1], 3, 1)
basis= [[1, 1, 0] [1, 0, 0]]
norm_rays = reshape([0, 0, 1], 3, 1)
test_normalise_rays(rays, basis, norm_rays)

