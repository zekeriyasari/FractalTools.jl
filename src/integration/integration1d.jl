# This file includes one dimensional fractal integration functions.

export integrate

# ------------------- Fractal Integration -------------------------------- #

"""
    $SIGNATURES

Returns the coefficients of the transformations of the ifs of the fractal interpolation.
"""
function coefficients(x, y, d)
    ifs = getifs(x, y, d) 
    a = [w.A[1, 1] for w in ifs.ws]
    c = [w.A[2, 1] for w in ifs.ws]
    e = [w.b[1] for w in ifs.ws]
    f = [w.b[2] for w in ifs.ws]
    (a, c, e, f)
end

function coefficients(x, y, z, d, h, l, m)
    ifs = gethiddenifs(x, y, z, d, h, l, m) 
    a = [w.A[1, 1] for w in ifs.ws]
    c = [w.A[2, 1] for w in ifs.ws]
    k = [w.A[3, 1] for w in ifs.ws]
    e = [w.b[1] for w in ifs.ws]
    f = [w.b[2] for w in ifs.ws]
    g = [w.b[3] for w in ifs.ws]
    (a, c, k, e, f, g)
end

"""
    $SIGNATURES

I = ∫f(x)dx

Numerically integrates the interpolation function for interpolation points (`x`, `y`). `x`
is a domain point, `y` is the range point. `d` is the vertical scaling factor (0., 1.)
"""
function integrate(x::AbstractVector, y::AbstractVector; d::Union{<:Real, <:AbstractVector}=0.001) 
    if d isa Real 
        d = fill(d, length(x) - 1)
    end 

    # Compute coefficients
    n = length(x) - 1
    a, c, e, f = coefficients(x, y, d)

    # Compute the integral
    α = sum(d.*a)
    K = zeros(n)
    for i=1:n
        g(x) = c[i] * x + f[i]
        K[i] = hquadrature(g, x[1], x[end]; rtol=1e-8, atol=0, maxevals=0)[1]
    end
    β = sum(a .* K)
    I = β / (1 - α)
    return I
end

"""
    integrate(x::AbstractVector{T}, y::AbstractVector{T}, g::Function;
        d=0.001) where {T<:Real}

Computes the numerical integration ``I = ∫g(u)f(u)du`` where `f` is the interpolation
function for the points (`x`, `y`). `g` is any single variable function. `d` is the vertical
scaling factor.
"""
function integrate(x::AbstractVector, y::AbstractVector, g::Function; d::Union{<:Real,<:AbstractVector}=0.001)
    # Size check 
    if d isa Real 
        d = fill(d, length(t) - 1) 
    end 

    # Compute coefficients
    n = length(x) - 1
    delta_x = diff(x)
    delta_y = diff(y)
    a = delta_x / (x[end] - x[1])
    c = delta_y / (x[end] - x[1]) .- d * (y[end] - y[1]) / (x[end] - x[1])
    e = x[end] / (x[end] - x[1]) * x[1 : end - 1] -
        x[1] / (x[end] - x[1]) * x[2 : end]
    f = x[end] / (x[end] - x[1]) * y[1 : end - 1] -
        x[1] / (x[end] - x[1]) * y[2 : end] .-
        d * (x[end] * y[1] - x[1] * y[end]) / (x[end] - x[1])

   K = zeros(n)
   for i = 1 : n
       g̃ = x̃ -> g(a[i] * x̃ + e[i]) * (c[i] * x̃ + f[i])
       K[i] = hquadrature(g̃, x[1], x[end]; rtol=1e-8, atol=0, maxevals=0)[1]
   end
   I = sum(abs.(a) .* K ./ (1 - d))
   return I
end

# ----------------------- Hidden Fractal Integration ------------------------- #

"""
    $SIGNATURES

I = ∫f(x)dx

Numerically integrates the interpolation function for interpolation points (`x`, `y`). `x`
is a domain point, `y` is the range point. `d` is the vertical scaling factor (0., 1.)
"""
function integrate(x::AbstractVector, y::AbstractVector, z::AbstractVector; 
    d::Union{<:Real, <:AbstractVector}=0.001,
    h::Union{<:Real, <:AbstractVector}=0.001,
    l::Union{<:Real, <:AbstractVector}=0.001,
    m::Union{<:Real, <:AbstractVector}=0.001) 

    d = isa(d, Real) ? fill(d, length(x) - 1) : d 
    h = isa(h, Real) ? fill(h, length(x) - 1) : h 
    l = isa(l, Real) ? fill(l, length(x) - 1) : l 
    m = isa(m, Real) ? fill(m, length(x) - 1) : m 

    # Compute coefficients
    n = length(x) - 1
    a, c, k, e, f, g = coefficients(x, y, z, d, h, l, m)

    # Compute the integral
    x0, xN = x[1], x[end] 
    Φ = [1 - sum(a .* d)     -sum(a .* h); 
        -sum(a .* l)        1 - sum(a .* m)] 
    θ = [sum(a .* ((c / 2 * (xN^2 - x0^2) + f * (xN - x0)))), 
        sum(a .* ((k / 2 * (xN^2 - x0^2) + g * (xN - x0))))]
    I = Φ \ θ
    return I[1]
end
