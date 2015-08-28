using Distributions

A = [[3 -4 1 0],
     [2 0 4 -1],
     [-4 -7 2 4],
     [-1 0 20 2],
     [6 -5 -4 2]]
B = [A, eye(Int, size(A,2))]

## C = EllipticalScenGen.remove_redundant_constraints(B)
## chernikova_general(C)
## chernikova_general(D)

## m, n = size(B)
## S = [[eye(eltype(B), n) B'], [-eye(eltype(B), n) -B']]
## ChernMat = EllipticalScenGen.ChernMat
## mat = ChernMat(S, m, n, 2*n)
## chernikova(mat, 1)

## using EllipticalScenGen
## q = ones(6)
## K = quota_cone(q)
## A = int(K.A)
## chernikova_general(A, 2)

# Function which checks anything generated from Chernikova output is in Polyhedral cone
function test_polyhedral_contains_finite_cone(A::Matrix{Int}, num_points::Int)
    X = chernikova(A)
    dim, gens = size(X)
    dist = Uniform(0.0, 100.0)
    B = [A; eye(Int,dim)]
    for i in 1:num_points
        λ = rand(dist, gens)
        @test all(B*(X*λ) .>= 0.0)
    end
end

# Function which checks that the polyhedral cone (intersected with the positive quadrant)
# is contained in the finitely generated cone output by the Chernikova function.
# This test is done by checking that the dual of the latter is contained in
# the dual of the former.
function test_finite_contains_polyhedral_cone(A::Matrix{Int}, num_points::Int)
    dim = size(A, 2)
    B = [A; eye(Int, dim)]  # Must add positivity constaints to original polyhedral cone
    dual_gens= size(B,1)
    X = chernikova(B)
    dist = Uniform(0.0, 100.0)
    for i in 1:num_points
        λ = rand(dist, dual_gens)
        @test all(X'B'λ .>= 0)
    end
end
    
test_polyhedral_contains_finite_cone(A, 10000)
test_finite_contains_polyhedral_cone(A, 10000)
test_polyhedral_contains_finite_cone(B, 10000)
test_finite_contains_polyhedral_cone(B, 10000)
