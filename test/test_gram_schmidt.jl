function test_gram_schmidt{T<:Integer}(A::Matrix{T})
    B = cones.gram_schmidt(A)
    @test cones.check_column_span_same(A, B)
    @test cones.check_columns_orthogonal(B)
end

srand(1)
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
