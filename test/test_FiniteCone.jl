# Tests whether points project correctly
using Test
using Cones

@testset "FiniteCone projection - #generators == dim == 2..." begin
    cone_1 = FiniteCone([[0.0, 1.0] [1.0, 1.0]])
    test_data_1 = [([-1.0, 3.0], [0.0, 3.0]),
                ([-1.0, -1.0], [0.0, 0.0]),
                ([2.0, 0.0], [1.0, 1.0]),
                ([0.25, 1.0], [0.25, 1.0])]
    for (x, proj_x) in test_data_1
        @test project(cone_1, x) ≈ proj_x
    end
end


@testset "FiniteCone projection - 3 = #generators > dim = 2..." begin
    cone_2 = FiniteCone([[0.0, 1.0] [1.0, 1.0] [0.5, 1.0]])
    test_data_2 = [([-1.0, 3.0], [0.0, 3.0]),
                ([-1.0, -1.0], [0.0, 0.0]),
                ([2.0, 0.0], [1.0, 1.0]),
                ([0.25, 1.0], [0.25, 1.0])]
    for (x, proj_x) in test_data_2
        @test project(cone_2, x) ≈ proj_x
    end
end

@testset "FiniteCone projection - 2 = #generators < dim = 3..." begin
    cone_3 = FiniteCone([[0.0, 0.0, 1.0] [1.0, 1.0, 1.0]])
    test_data_3 = [([-1.0, -1.0, 1.0], [0.0, 0.0, 1.0]),
                ([3.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
                ([-1.0, -1.0, -1.0], [0.0, 0.0, 0.0]),
                ([1.0, 2.0, 3.0], [1.5, 1.5, 3.0]),
                ([2.5, 2.5, 6.0], [2.5, 2.5, 6.0])]
    for (x, proj_x) in test_data_3
        @test project(cone_3, x) ≈ proj_x
    end
end

@testset "FiniteCone membership # generations = 3, dim = 3" begin
    cone_4 = FiniteCone([[1.0, 1.0, 0.0] [1.0, 0.0, 1.0] [0.0, 1.0, 1.0]])
    test_data_4 = [([1.0, 1.0, 0.0], true),
                    ([1.0, 1.0, 1.0], true),
                    ([1.0, 0.0, 0.0], false)]
    for (x, res) in test_data_4
        @test (x ∈ cone_4) == res
    end
end

@testset "FiniteCone - dual - 2D" begin
    K1 = FiniteCone([[1, 2] [2, 1]]) # Polyhedral form is (y <= 2x, y >= (1/2)x )
    K2 = Cones.dual(K1) # (y >= -2x, 2y >= -x)
    pts = [[1,1], [-1,2], [2,-1]] # Points in dual cone
    for i in 1:10
        c = K1.A*rand(2) # random point in original cone
        @test c'pts[1] >= 0
        @test c'pts[2] >= 0
        @test c'pts[3] >= 0
    end
end