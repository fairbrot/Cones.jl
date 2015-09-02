import Base: num, den

num{T<:Rational}(x::Array{T}) = map(y->y.num, x)
den{T<:Rational}(x::Array{T}) = map(y->y.den, x)

# Transform a rational matrix into an integer matrix by multiplying each row
intmat{T<:Rational}(A::Matrix{T}) =  mapslices(z->int(z*(lcm(den(z))//gcd(num(z)))), A, 2)

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
        if (g = gcd(A[:,i])) > 1; norm_A[:, i] = broadcast(div, norm_A[:, i], g); end
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
