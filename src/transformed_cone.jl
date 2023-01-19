struct TransformedCone <: Cone
    cone::Cone # Original cone
    U::Matrix{Float64} # Transformation
    function TransformedCone(cone::Cone, U::Matrix{Float64})
        n = length(cone)
        @assert size(U,1) == size(U,2) "Matrix must be square"
        @assert size(U,1) == n "Matrix and cone must have compatible dimensions"
        @assert maximum(abs.(Matrix(1.0I, n, n) - U'U)) < 1e-6 "Matrix must be orthonormal"
        return new(cone, U)
    end
end

∈(p::Vector{T}, K::TransformedCone) where {T<:Real} = K.U'p ∈ K.cone

dual(K::TransformedCone) = TransformedCone(dual(K.cone), K.U)