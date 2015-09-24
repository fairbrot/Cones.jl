# Checks output from chernikova_general function an input matrix
# defines a specified cone. That is,
# { x: Ax ≥ 0 } == { B y + C z : y ≥ 0, z free }
#
# Arguments:
# A::Matrix - constraint matrix (rows correspond to individual constraints)
# bidrays::Matrix - matrix of unidirectional rays (columns correspond to individual rays)
# unirays::Matrix - matrix of (normalised) bidirectional rays (columns correspond to individual rays)
function test_chernikova_general{T1<:Union(Integer, Rational), T2<:Integer}(A::Matrix{T1}, bidrays::Matrix{T2}, rays::Matrix{T2})
    norm_rays = cones.normalise_rays(rays, bidrays)
    chern_bidrays, chern_rays = chernikova_general(A)
    norm_chern_rays = cones.normalise_rays(chern_rays, chern_bidrays)
    @test cones.check_columns_same(norm_rays, norm_chern_rays)
    @test cones.check_column_span_same(bidrays, chern_bidrays)
end

# Pointed 2D cone
println("\tExample 1...")
A = [[-1//5 4//5],
     [4//5 -1//5]]
rays = [[1,4] [4, 1]]
bidrays = Array(Int, 2, 0)
test_chernikova_general(A,bidrays, rays)

# Hyperplane
println("\tExample 2...")
A = [[-1 -1 0], [1 1 0]]
rays = Array(Int, 3, 0)
bidrays = [[1, -1, 0] [0, 0, 1]]
test_chernikova_general(A, bidrays, rays)

# 1D lineality space
println("\tExample 3...")
A = [[-3 3 0], [2 2 0]]
rays = [[-1, 1, 1] [2, 2, -1]]
bidrays = [0 0 1]'
test_chernikova_general(A, bidrays, rays)

## function unirays{T<:Integer}(bidrays::Matrix{T}, unirays::Matrix{T})
##     d, n1 = size(bidrays)
##     n2 = size(unirays, 2)
##     rays = Array(T, d, 2*n1 + n2)
##     for i in 1:n1
##         rays[:, 2*i - 1] = bidrays[:,i]
##         rays[:, 2*i] = -bidrays[:,i]
##     end
##     rays[:, (2*n1 + 1):(2*n1 + n2)] = unirays
##     return rays
## end

## # Function which checks anything generated from Chernikova output is in Polyhedral cone
## function test_polyhedral_contains_finite_cone(A::Matrix{Int}, num_points::Int)
##     X = chernikova(A)
##     dim, gens = size(X)
##     dist = Uniform(0.0, 100.0)
##     B = [A; eye(Int,dim)]
##     for i in 1:num_points
##         λ = rand(dist, gens)
##         @test all(B*(X*λ) .>= 0.0)
##     end
## end
