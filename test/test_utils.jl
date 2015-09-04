function test_gram_schmidt{T<:Integer}(A::Matrix{T})
    B = cones.gram_schmidt(A)
    @test cones.check_column_span_same(A, B)
    @test cones.check_columns_orthogonal(B)
end

srand(1)

println("\t...gram_schmidt")

N = 10         # Trials
m, n = 5, 3
counter = 0
while counter < N
    # Randomly generate full column rank matrix
    A = rand(Int, m, n) % 5
    if rank(A) < n continue end
    test_gram_schmidt(A)
    counter += 1
end

# Test normalise_rays
# 
function test_normalise_rays{T<:Integer}(rays::Matrix{T}, basis::Matrix{T}, norm_rays::Matrix{T})
    orth_basis = cones.gram_schmidt(basis)
    norm_rays_out = cones.normalise_rays(rays, basis)
    @test cones.check_columns_same(norm_rays, norm_rays_out)
end

println("\t...normalise_rays")

# 0D lineality space
rays = [[2, 3, 1] [-1, 3, 2]]
basis = Array(Int, 3, 0)
norm_rays = rays
test_normalise_rays(rays, basis, norm_rays)

# 1D lineality space
rays = [[2, 0, 2] [0, 1, 1]]
basis= [0 0 1]'
norm_rays = [[1, 0, 0] [0, 1, 0]]
test_normalise_rays(rays, basis, norm_rays)

rays = [[0, 2, 1] [3, 1, 1]]
basis = [-3 -3 0]'
norm_rays = [[-1, 1, 1] [1, -1, 1]]
test_normalise_rays(rays, basis, norm_rays)
         
# 2D lineality space
rays = [1 1 1]'
basis= [[1, 1, 0] [1, 0, 0]]
norm_rays = [0 0 1]'
test_normalise_rays(rays, basis, norm_rays)

