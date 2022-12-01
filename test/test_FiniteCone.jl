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