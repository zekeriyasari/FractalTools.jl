# This file includes hidden ifs 

export hiddenfif, gethiddenifs, wraphidden

"""
    $(SIGNATURES)

Interpolates the data pairs (xi, yi) for xi ∈ `x` and yi ∈ y. `f0` is the initial function and `niter` is the number of iterations.
"""
function hiddenfif(x::AbstractVector, y::AbstractVector, z::AbstractVector;  
    d::Union{<:Real,<:AbstractVector}=0.01, h::Union{<:Real,<:AbstractVector}=0.01, 
    l::Union{<:Real,<:AbstractVector}=0.01, m::Union{<:Real,<:AbstractVector}=0.01, 
    f0=x->[0., 0.], niter::Int=5)

    # Check types.
    d = isa(d, Real) ? fill(d, length(x) - 1) : d 
    h = isa(h, Real) ? fill(h, length(x) - 1) : h 
    l = isa(l, Real) ? fill(l, length(x) - 1) : l 
    m = isa(m, Real) ? fill(m, length(x) - 1) : m 

    # Construct IFS 
    ifs = gethiddenifs(x, y, z, h, d, l, m)

    # Construct domains 
    domains = getintervals(x) 

    # Construct itp 
    itp = ∘((wraphidden for i in 1 : niter)...)

    # Return interpolant 
    Interpolant(ifs, domains, itp((f0, ifs, domains))[1])
end

"""

    $(SIGNATURES)

Returns the IFS corresponding to interpolation data `x`, `y`. `d` is the scaling parameters. See the notes on [docs](https://zekeriyasari.github.io/FractalTools.jl/stable/manual/one_dimensional_interpolation/#D-Hidden-Fractal-Interpolation)
"""
function gethiddenifs(x::AbstractVector, y::AbstractVector, z::AbstractVector,  
    d::AbstractVector, h::AbstractVector, l::AbstractVector, m::AbstractVector)
    K = 1 / (x[1] - x[end])
    Ω = [1 -1; -x[end] x[1]]
    X = [x[1 : end - 1, :]'; x[2 : end, :]']
    Y = [y[1 : end - 1, :]'; y[2 : end, :]'] - [y[1] * d'; y[end] * d'] -[z[1] * h'; z[end] * h'] 
    Z = [z[1 : end - 1, :]'; z[2 : end, :]'] - [y[1] * l'; y[end] * d'] -[z[1] * m'; z[end] * h'] 
    res = K * map(ω -> Ω * ω, hcat.(eachcol(X), eachcol(Y), eachcol(Z)))
    ws = map(enumerate(res)) do item 
        n, (a, e, c, f, k, g) = item 
        Transformation([a 0 0 ; c d[n] h[n]; k l[n] m[n]], [e, f, g])
    end
    IFS(ws)
end

"""
   $(SIGNATURES)
   
Functional that takes initiali function `func` and returns a wrapped function.
"""
function wraphidden((func, ifs, domains))
    function wrapped(x)
        n = findfirst(I -> x ∈ I, domains)
        w = ifs.ws[n]
        (a, c, k, _, d, l, _, h, m), (e, f, g) = w.A, w.b
        ω = (x - e) / a 
        [d h; l m] * func(ω) + [c; k] * ω + [f, g]
    end, ifs, domains
end
