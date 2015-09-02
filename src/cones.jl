module cones

VERSION < v"0.4-" && using Docile

using MathProgBase, Gurobi

import Base: length

export chernikova, chernikova_general, Cone, FiniteCone, PolyhedralCone, project
# package code goes here
include("utils.jl")
include("LCP_julia.jl")
include("chernikova.jl")
include("chernikova_general.jl")
include("types.jl")

end # module
