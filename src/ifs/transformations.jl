# This file includes the tranformations.

export Transformation
export dimension, contfactor

""" 
    $(TYPEDEF) 

Affine transformation 

# Fields 

    $(TYPEDFIELDS)
"""
struct Transformation{T1, T2}
    "Matrix of the transformation"
    A::T1 
    "Vector of the transformation"
    b::T2 
    function Transformation(A::AbstractMatrix, b::AbstractVector)
        size(A, 2) == length(b) || throw(DimensionMismatch("Size of `A` is not compatible with `b`"))
        new{typeof(A), typeof(b)}(A, b)
    end 
end 

(w::Transformation)(x) = w.A * x + w.b

"""
    $SIGNATURES 

Returns dimension of `w`.
"""
dimension(w::Transformation) = size(w.A,1)

"""
    $SIGNATURES 

Returns contraction factor of `w`. Contraction factor is computed as the norm of `w.A`.
"""
contfactor(w::Transformation) = norm(w.A)



