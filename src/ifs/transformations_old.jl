# This file contains Transforms of FractalTools

using LinearAlgebra

"""
    AbstractAffineTransform <: Function

Abstract type for affine transformations. An `AbstractAffineTransform` transforms a vector `x` to ` A * x + b` where A
is a matrix and b is a vector.
"""
abstract type AbstractAffineTransform <: Function end

"""
    AffineTransform(A, b)

Constructs and affine transform with `A` and `b`.
"""
struct AffineTransform <: AbstractAffineTransform
    A::Matrix{Float64}
    b::Vector{Float64}
    function AffineTransform(A, b)
        size(A, 2) == length(b) ? new(A, b) : error("Matrix and vector dimension mismatch.")
    end
end

"""
    AffineContraction(A, b, contfactor)

Constructs an `AffineContraction` with `A` and `b`. `contfactor` is the contration factor of `AffineContraction`. The contraction is defined as
```math
    ‖ f(x) - f(y) ‖ ≤ α ‖ x- y ‖  ∀x,y ∈ Χ
```
where ``f`` is and `AffineContraction` and ``α`` is the `contfactor`.
"""
struct AffineContraction <: AbstractAffineTransform
    A::Matrix{Float64}
    b::Vector{Float64}
    contfactor::Float64
end
AffineContraction(A, b) = AffineContraction(A, b, norm(A))

show(io::IO, trfm::AffineTransform) = print(io, "AffineTransform(A:$(trfm.A), b:$(trfm.b))")
show(io::IO, trfm::AffineContraction) = print(io, "AffineContraction(A:$(trfm.A), b:$(trfm.b), contfactor:$(trfm.contfactor))")
display(trfm::AffineTransform) = println("AffineTransform(A:$(trfm.A), b:$(trfm.b))")
display(trfm::AffineContraction) = println("AffineContraction(A:$(trfm.A), b:$(trfm.b), contfactor:$(trfm.contfactor))")

(trfm::AbstractAffineTransform)(x) = trfm.A * x + trfm.b

"""
    dimension(trfm::AbstractAffineTransform)

Returns the dimension of `trfm`.
"""
dimension(trfm::AbstractAffineTransform) = size(trfm.A, 1)

"""
    contfactor(trmf::AffineTransform)

Returns contraction factor of `trfm`. Conracttion factor is computed as the norm of `trfm.A`.
"""
contfactor(trmf::AffineTransform) = norm(trmf.A)

"""
    contfactor(trfm::AffineContraction)

Returns the contraction factor of `trfm`.
"""
contfactor(trfm::AffineContraction) = trfm.contfactor
