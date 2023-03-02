import Base: ∈
using LinearAlgebra: norm, dot


"""
# Description
Type representing a scaled second order cone:
{ x ∈ R^d : x_1 ≥ α ||x_{2:n} || }

# Arguments
`d::Int`: dimension of cone
`α::Float64`: scaling parameter
"""
struct SOCone <: Cone
    d::Int # Dimension of space
    α::Float64 # Scaling parameter
end

function ∈(p::Vector{T}, cone::SOCone) where {T<:Real}
    return p[1] >= cone.α*norm(view(p, 2:length(p)))
end

length(cone::SOCone) = cone.d

dual(cone::SOCone) = SOCone(cone.d, 1/cone.α)

"Checks whether `p` is on the boundary of `cone`. Provided for testing purposes."
function on_boundary(p::Vector{T}, cone::SOCone) where {T<:Real}
    return p[1] ≈ cone.α*norm(view(p, 2:length(p)))
end

"Returns a point `c` such that c[i] = p[i] and `c` is on boundary of `cone`. Provided for testing purposes."
function move_to_boundary(p::Vector{T}, cone::SOCone) where {T<:Real}
    c = copy(p)
    c[1] = cone.α*norm(view(p, 2:length(p)))
    return c
end

struct GeneralSOCone
    u::Vector{Float64} # cone axis
    α::Float64
    function GeneralSOCone(u::Vector{T}, α::Real) where {T<:Real} 
        @assert norm(u) ≈ 1.0 
        return new(u, α)
    end
end

length(cone::GeneralSOCone) = length(cone.α)


function ∈(p::Vector{T}, cone::GeneralSOCone) where {T<:Real}
    a = dot(p, cone.u)
    proj = p - a*cone.u
    return a+sqrt(eps()) >= cone.α*norm(proj)
end

dual(cone::GeneralSOCone) = GeneralSOCone(cone.u, 1/cone.α)

