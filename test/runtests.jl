using cones
using Base.Test

println("Running test_FiniteCone.jl...")
include("test_FiniteCone.jl")

println("Running test_chernikova.jl...")
include("test_chernikova.jl")

println("Running test_PolyhedralCone.jl...")
include("test_PolyhedralCone.jl")
