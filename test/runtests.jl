using Cones
using Test

println("Running test_utils.jl...")
include("test_utils.jl")

println("Running test_FiniteCone.jl...")
include("test_FiniteCone.jl")

println("Running test_chernikova.jl...")
include("test_chernikova.jl")

println("Running test_chernikova_general.jl...")
include("test_chernikova_general.jl")

println("Running test_PolyhedralCone.jl...")
include("test_PolyhedralCone.jl")

println("Running test_polytope.jl...")
include("test_polytope.jl")
