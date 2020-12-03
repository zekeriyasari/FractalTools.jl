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
    A::T1 
    b::T2 
    Transformation(A::AbstractMatrix, b::AbstractVector) = new{typeof(A), typeof(b)}(A, b)
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



