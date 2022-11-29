import Base.length, Base.∈

abstract type Cone end

"""
# Description
Type representing a finitely generated cone. This is the positive
hull of a finite collection of vectors a₁, … , aₙ:

{ ∑ λᵢ xᵢ : λ ≥ 0 }

# Arguments
`A::Matrix{Float64}`: Matrix of cone generators, where each column corresponds to a generator
"""
struct FiniteCone <: Cone
    A::Matrix{Float64}    # Matrix where each column is a cone generator
    AtA::Matrix{Float64}  # Cross-product of cone generator matrix
    num_gen::Int64
    FiniteCone(A::Matrix{T}) where T<:Real = new(float(A), float(A'A), size(A,2))
end

function project(cone::FiniteCone, p::Vector{Float64})
    w,z = lcp_solve(cone.AtA, -cone.A'p)
    cone.A * z
end

function ∈(p::Vector{Float64}, cone::FiniteCone; optimizer=Gurobi.Optimizer)
    n = length(p)
    m = cone.num_gen
    model = Model(optimizer)
    set_silent(model)
    @variable(model, λ[1:m] ≥ 0)
    @constraint(model, [i=1:n], sum(cone.A[i,j]*λ[j] for j in 1:m) == p[i])
    optimize!(model)
    stat = termination_status(model)
    if stat == MOI.OPTIMAL
        return true
    elseif stat == MOI.INFEASIBLE
        return false
    else
        throw(ErrorException("Optimization status should be optimal or infeasible"))
    end
end

# Dimension of ambient space
length(cone::FiniteCone) = size(cone.A, 1)

"""
# Description
Type representing a Polyhedra cone. This is the
intersection of a finite collection of half-spaces

{x : aᵢ x ≥ 0 ∀ i} = { x: Ax ≥ 0 }

# Arguments
`A::Matrix{Float64}`: Constraint matrix of polyhedral cone
"""
struct PolyhedralCone{T<:Real} <: Cone
    A::Matrix{T}
    n::Int   # Dimension
    m::Int   # Number of constraints
    function PolyhedralCone{T}(A::Matrix{T}) where T<: Real
        return new{T}(A, size(A,2), size(A,1))
    end
end

function PolyhedralCone(A::Matrix{T}, check=false) where T<:Integer
    if check
        A = remove_redundant_constraints(A)
    end
    PolyhedralCone{T}(A)
end

function project(cone::PolyhedralCone, p::Vector{Float64}; optimizer=Gurobi.Optimizer) 
    n = length(p)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:n])
    @constraint(model, cone.A * x .>= 0)
    @objective(model, Min, sum(x[i]*x[i] - 2*p[i]*x[i] for i in 1:n))
    optimize!(model)
    @assert termination_status(model) == MOI.OPTIMAL
    return value.(x)
end

# Dimension of ambient space
length(cone::PolyhedralCone) = size(cone.A, 2)
