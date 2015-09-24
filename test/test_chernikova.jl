using Distributions

# Checks output from chernikova function for some input matrix
# defines a specified finitely generated cone. That is,
# { x: Ax ≥ 0 } == { B y : y ≥ 0 }
#
# Arguments:
# A::Matrix - constraint matrix (rows correspond to individual constraints)
# B::Matrix - matrix of cones rays (columns correspond to individual rays)
function test_chernikova{T1<:Union(Integer, Rational), T2<:Integer}(A::Matrix{T1}, B::Matrix{T2})
    rays = cones.norm_cols(B)
    chern_rays = chernikova(A)
    @test cones.check_columns_same(rays, chern_rays)
end

# Simple two rays example
println("\tExample 1...")
A = [[-1//5 4//5],
     [4//5 -1//5]]
B = [[1,4] [4, 1]]
test_chernikova(A,B)


# Redundant constraints
println("\tExample 2...")
A = [[1 1 0],
     [1 0 1],
     [1 1 1],
     [1 1 1],
     [0 1 1]]
B = eye(Int, 3)
test_chernikova(A,B)

# Empty set example
println("\tExample 3...")
A = [[-1 -1 -1]]
B = Array(Int, 3, 0)
test_chernikova(A, B)


# One ray example
println("\tExample 4...")
A = [[1 -1], [-1 1]]
B = [1 1]'
test_chernikova(A, B)


# Example from paper
# "Algorithm for finding a general formula for the non-negative solutions of a system of linear inequalities"
# by Chernikova (1964)
println("\tExample 5...")
A = [[3 -4 1 0],
     [2 0 4 -1],
     [-4 -7 2 4],
     [-1 0 20 2],
     [6 -5 -4 2]]

B = [[1, 0, 0, 2] [0, 0, 1, 4] [0, 1, 4, 16] [0, 0, 2, 4]  [0, 32, 128, 336] [6, 0, 0, 6] [21, 12, 0, 42] [52, 40, 4, 120] [20, 0, 32, 4] [24080, 25760, 30800, 53760]]
test_chernikova(A, B)
