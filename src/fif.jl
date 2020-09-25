# This file includes the tools for one - dimensional fractal interpolation 

export Interval, AbstractInterpolant, Interpolant, fif, getifs, getintervals
import Base.∈

abstract type AbstractInterpolant end

"""
    $(TYPEDEF)

Real closed interval

# Fields 

    $(TYPEDFIELDS)
"""
struct Interval{T<:Real}
    lb::T 
    ub::T 
end 
Interval(lb::Real, ub::Real) = Interval(promote(lb, ub)...)

"""
    $(SIGNATURES)

Returns true if `x` is in `I`.

# Example 
```jldoctest 
julia> I = Interval(1, 2)
Interval{Int64}(1, 2)

julia> 1.5 ∈ I 
true
```
"""
∈(x::Real, I::Interval) = I.lb ≤ x ≤ I.ub


"""
    $(TYPEDEF)

One dimensional fractal interpolant 

# Fields 

    $(TYPEDFIELDS)
"""
struct Interpolant{T1<:IFS, T2<:AbstractVector{<:Interval}, T3} <: AbstractInterpolant
    ifs::T1         # IFS of interpolations. Transformations are computed from the data 
    intervals::T2   # Just-touching intervals 
    itp::T3         # Interpolant function 
end 

(interp::Interpolant)(x) = interp.itp(x)

"""
    $(SIGNATURES)

Interpolates the data pairs (xi, yi) for xi ∈ `x` and yi ∈ y. `f0` is the initial function and `niter` is the number of iterations.
"""
function fif(x::AbstractVector, y::AbstractVector, d::AbstractVector; f0=zero, niter::Int=5)
    # Construct IFS 
    ifs = getifs(x, y, d)

    # Construct intervals 
    intervals = getintervals(x) 

    # Construct itp 
    itp = ∘((wrap for i in 1 : niter)...)

    # Return interpolant 
    Interpolant(ifs, intervals, itp((f0, ifs, intervals))[1])
end

"""

    $(SIGNATURES)

Returns the IFS corresponding to interpolation data `x`, `y`. `d` is the scaling parameters.
"""
function getifs(x::AbstractVector, y::AbstractVector, d::AbstractVector)
    K = 1 / (x[1] - x[end])
    Ω = [1 -1; -x[end] x[1]]
    X = [x[1 : end - 1, :]'; x[2 : end, :]']
    Y = [y[1 : end - 1, :]'; y[2 : end, :]'] - [y[1] * d'; y[end] * d']
    res = K * map(ω -> Ω * ω, hcat.(eachcol(X), eachcol(Y)))
    ws = map(enumerate(res)) do item 
        n, (a, e, c, f) = item 
        Transformation([a 0; c d[n]], [e, f])
    end
    IFS(ws)
end

"""
    $(SIGNATURES)

Returns the just-touching intervals 
"""
getintervals(x::AbstractVector) = map(n -> Interval(x[n], x[n + 1]), 1 : length(x) - 1)


"""
   $(SIGNATURES)
   
Functional that takes initiali function `func` and returns a wrapped function.
"""
function wrap((func, ifs, intervals))
    function wrapped(x)
        n = findfirst(I -> x ∈ I, intervals)
        w = ifs.ws[n]
        (a, c, _, d), (e, f) = w.A, w.b
        ω = (x - e) / a 
        c * ω + d * func(ω) + f
    end, ifs, intervals
end
