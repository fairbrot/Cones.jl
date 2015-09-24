# Find generators for the polytope
# P = { x : Ax ≥ b }
# that is, calculate the following representation
# P = { ∑ᵢ λᵢUᵢ + ∑ⱼ μⱼVⱼ + ∑ₖ νₖWₖ : ∑ᵢλᵢ = 1, μ ≥ 0 }
#
# The matrices U, V and W have columns which
# represent the extremal points, unidirectional rays and bidirectional rays
# of the polytope respectively.
#
# Returns:
# (points, rays, bidrays):
#   `points::Matrix{Rational{T}}`
#   `rays::Matrix{T}`
#   `bidrays::Matrix{T}`
function find_generators{T<:Integer}(A::Matrix{T}, b::Vector{T})
    # Create homogenized system
    m, n = size(A)
    Ab = [[A -b], zeros(T, 1, n+1)]
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
            rays[:,r_count] = uni[1:n,j]//uni[n+1,j]
            r_count+=1
        end
    end

    bidrays = bid[1:n,:]
    return points, rays, bidrays
end
