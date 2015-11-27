type Polytope{T<:Integer}
    points::Matrix{Rational{T}}
    rays::Matrix{T}
    bidrays::Matrix{T}
end

function find_generators{T<:Rational}(A::Matrix{T}, b::Vector{T})
    mat = [A b]
    imat = intmat(mat)
    Aint, bint = imat[:,1:end-1], imat[:,end]
    find_generators(Aint, bint)
end

@doc """
# Description
Find generators for the polytope
  P = { x : Ax ≥ b }
that is, calculate the following representation
  P = { ∑ᵢ λᵢUᵢ + ∑ⱼ μⱼVⱼ + ∑ₖ νₖWₖ : ∑ᵢλᵢ = 1, μ ≥ 0 }

This is done through the use of the generalised Chernikova
algorithm on the cone homogenization of the above system.

# Returns
(points, rays, bidrays):
  `points::Matrix{Rational{T}}`: matrix whose columns correspond extremal points of polytope
  `rays::Matrix{T}`: matrix whose columns correspond to unidirectional rays of polytope
  `bidrays::Matrix{T}`: matrix whose columns correspond to bidirectional rays of polytopes
""" ->
function find_generators{T<:Integer}(A::Matrix{T}, b::Vector{T})
    # Create homogenized system
    m, n = size(A)
    Ab = [[A -b]; zeros(T, 1, n+1)]
    Ab[m+1,n+1] = one(T)
    bid, uni = chernikova_general(Ab)

    n_bid = size(bid, 2)
    n_uni = size(uni, 2)
    
    # Count number of trailing zeros (rays of polytope)
    # Note that all bidirectional rays in homogenized system
    # must have a trailing zero
    # (due to positivity constraint on this coordinate)
    n_zeros_uni = 0
    for j in 1:n_uni
        if uni[n+1,j] == 0 n_zeros_uni += 1 end
    end

    # If no unidirectional rays intersect with {x : x[n+1] == 0}
    # then polytope is empty
    if n_zeros_uni == n_uni
        return Array(Rational{T}, n, 0), Array(T, n, 0), Array(T, n, 0)
    end
    
    # Allocate matrices and get points.
    # The vertices of the polytope correspond to the unidirectional rays with non-zero trailing coordinate
    points = Array(Rational{T}, n, (n_uni - n_zeros_uni))
    rays = Array(T, n, n_zeros_uni)
    p_count = 1
    r_count = 1
    for j in 1:n_uni
        if uni[n+1, j] != 0
            points[:, p_count] = uni[1:n,j]//uni[n+1,j]
            p_count+=1
        else
            rays[:,r_count] = uni[1:n,j]
            r_count+=1
        end
    end

    bidrays = bid[1:n,:]
    return points, rays, bidrays
end

@doc """
# Description
Calculates the minimum and maximum values of
the projection of the polytope onto
specified vectors:

  min_{x ∈ P} Tᵢᵗx 
  max_{x ∈ P} Tᵢᵗx

where Tᵢ are the rows of a matrix T.
# Arguments
* `P::Polytope{I}`: Polytope to project
* `T::Matrix{I}`: Matrix onto whose columns we project polytope
# Returns
`(mins::Array{Rational{T}}, maxs::Array{Rational{T}})`
""" ->
function min_max_projections{I<:Integer}(P::Polytope{I}, T::Matrix{Rational{I}} = eye(Rational{I}, length(P)))
    n = size(T, 2)
    n_points, n_rays, n_bidrays = size(P.points, 2), size(P.rays, 2), size(P.bidrays, 2)
    mins, maxs = fill(typemax(Rational{I}), n), fill(typemin(Rational{I}), n)
    for i in 1:n
        for j in 1:n_points
            tp = dot(T[:,i], P.points[:,j])
            if tp < mins[i]
                mins[i] = tp
            end
            if tp > maxs[i]
                maxs[i] = tp
            end
        end
        for j in 1:n_rays
            tr = dot(T[:,i], P.rays[:,j])
            if tr < 0
                mins[i] = -Inf
                break
            elseif tr > 0
                maxs[i] = Inf
                break
            end
        end
        for j in 1:n_bidrays
            if dot(T[:,i], P.bidrays[:,j]) != 0
                maxs[i] = Inf
                mins[i] = -Inf
                break
            end
        end
    end
    return mins, maxs
end

length{T<:Integer}(P::Polytope{T}) = size(P.points, 1)
