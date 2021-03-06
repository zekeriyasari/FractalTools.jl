# This file includes the tools for one - dimensional fractal interpolation 

export Interval, AbstractInterpolant, Interpolant, interpolate, getifs, getintervals
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

julia> 1.5 ∈ Interval(1, 2)
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
struct Interpolant{T1<:IFS, T2<:AbstractVector, T3} <: AbstractInterpolant
    ifs::T1         # IFS of interpolations. Transformations are computed from the data 
    domains::T2     # Just-touching domains 
    itp::T3         # Interpolant function 
end 

(interp::Interpolant)(args...) = interp.itp(args...)

"""
    $(SIGNATURES)

Interpolates the data pairs (xi, yi) for xi ∈ `x` and yi ∈ y. `f0` is the initial function and `niter` is the number of iterations.
"""
function interpolate(x::AbstractVector, y::AbstractVector; d::Union{<:Real, <:AbstractVector}, f0=zero, niter::Int=5)
    if d isa Real 
        d = fill(d, length(x) - 1)
    end 

    # Construct IFS 
    ifs = getifs(x, y, d)

    # Construct domains 
    domains = getintervals(x) 

    # Construct itp 
    itp = ∘((wrap for i in 1 : niter)...)

    # Return interpolant 
    Interpolant(ifs, domains, itp((f0, ifs, domains))[1])
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

Returns the just-touching domains 
"""
getintervals(x::AbstractVector) = map(n -> Interval(x[n], x[n + 1]), 1 : length(x) - 1)


"""
   $(SIGNATURES)
   
Functional that takes initiali function `func` and returns a wrapped function.
"""
function wrap((func, ifs, domains))
    function wrapped(x)
        n = findfirst(I -> x ∈ I, domains)
        w = ifs.ws[n]
        (a, c, _, d), (e, f) = w.A, w.b
        ω = (x - e) / a 
        c * ω + d * func(ω) + f
    end, ifs, domains
end
