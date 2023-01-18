module Cones

using LinearAlgebra: I, rank, dot
using JuMP
using Gurobi
using Printf: @printf

import Base: length


export chernikova, chernikova_general, Cone, FiniteCone, PolyhedralCone, project, Polytope, find_generators, min_max_projections, SOCone

# package code goes here
include("utils.jl")
include("LCP_julia.jl")
include("chernikova.jl")
include("chernikova_general.jl")
include("types.jl")
include("polytope.jl")
include("soc.jl")

end # module
