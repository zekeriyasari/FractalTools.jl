# This file includes the tranformations.

export Transformation

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