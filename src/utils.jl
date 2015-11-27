import Base: num, den

num{T<:Rational}(x::Array{T}) = map(y->y.num, x)
den{T<:Rational}(x::Array{T}) = map(y->y.den, x)

# Transform a rational matrix into an integer matrix by multiplying each row
intmat{T<:Rational}(A::Matrix{T}) =  mapslices(z->round(Int, z*(lcm(den(z))//gcd(num(z)))), A, 2)

function collinear(a::Array{Int}, b::Array{Int})
    length(a) == length(b) || error("Input arrays must have the same dimenisions")
    a_n, b_n = a/gcd(a), b/gcd(b)
    if a_n==b_n || a_n==-b_n; return true; end
    return false
end

# Normalises the columns of a matrix
function norm_cols(A::Matrix{Int})
    norm_A = copy(A)
    for i in 1:size(A,2)
        if (g = gcd(A[:,i])) > 1; norm_A[:, i] = div(norm_A[:, i], g); end
    end
    return norm_A
end

# Check matrices are the same up to permuting columns
function check_columns_same{T}(A::Matrix{T}, B::Matrix{T})
    m,n  = size(A)
    if size(B) != (m,n) return false end
    sA = Set{Array{T, 1}}()
    sB = Set{Array{T, 1}}()
    for i in 1:n
        push!(sA, A[:, i])
        push!(sB, B[:, i])
    end
    return sA == sB
end

# Check columns of two matrices span same linear subspace
# It is assumed columns in each matrix are linearly independent
function check_column_span_same{T<:Real}(A::Matrix{T}, B::Matrix{T})
    m, n = size(A)
    if size(B) != (m,n) return false end
    rankA = rank(A)
    return rank(A) == rank([A B])
end

# Tests whether all elements of vector are zero
function is_zero{T<:Real}(x::Vector{T})
    z = zero(T)
    for xi in x
        if xi != z return false end
    end
    return true
end

# Verify columns of input matrix are orthogonal (no zero columns allowed)
function check_columns_orthogonal{T<:Real}(A::Matrix{T})
    m, n = size(A)
    for i in 1:n
        if is_zero(A[:,i]) return false end
        for j in 1:(i-1)
            if dot(A[:,i], A[:,j]) != zero(T) return false end
        end
    end        
    return true
end

# Use a modified version of Gram-Schmitd to orthonalize columns of integer matrix
function gram_schmidt{T<:Integer}(A::Matrix{T})
    m, n = size(A)
    B = Array(T, m, n)
    for i in 1:n
        sq_norms = T[dot(B[:,j], B[:,j]) for j in 1:(i-1)]
        prod_sq_norms = prod(sq_norms)
        B[:,i] = prod_sq_norms * A[:,i]
        for j in 1:(i-1)
            B[:,i] -= div(prod_sq_norms,sq_norms[j]) * dot(A[:,i], B[:,j]) * B[:,j]
        end
        if (g = gcd(B[:,i])) > 1; B[:, i] = div(B[:, i], g); end
    end
    return B
end

# Transforms rays so that they are orthogonal to all elements in a linear subspace
# Arguments:
# `A::Matrix{T}`: matrix of rays (each column is a ray)
# `B::Matrix{T}`: matrix whose columns form a basis of linear subspace
function normalise_rays{T<:Integer}(A::Matrix{T}, B::Matrix{T})
    m, n = size(A)
    k = size(B,2) # Size of basis
    new_rays = Array(T, m, n)
    B = gram_schmidt(B) # must have orthogonal basis for linear subspace
    sq_norms = T[dot(B[:,j], B[:,j]) for j in 1:k]
    prod_sq_norms = prod(sq_norms)
    for i in 1:n
        new_rays[:,i] = prod_sq_norms * A[:,i]
        for j in 1:k
            new_rays[:,i] -= div(prod_sq_norms,sq_norms[j]) * dot(A[:,i], B[:,j]) * B[:,j]
            if (g = gcd(new_rays[:,i])) > 1; new_rays[:, i] = div(new_rays[:, i], g); end
        end
    end
    return new_rays
end
