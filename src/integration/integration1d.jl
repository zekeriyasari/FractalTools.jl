# This file includes one dimensional fractal integration functions.

"""
    integrate(x::AbstractVector, y::AbstractVector; d::Real=0.001)

I = ∫f(x)dx

Numerically integrates the interpolation function for interpolation points (`x`, `y`). `x`
is a domain point, `y` is the range point. `d` is the vertical scaling factor (0., 1.)
"""
function integrate(x::AbstractVector, y::AbstractVector; d::Real=0.001) 
    # Compute coefficients
    n = length(x) - 1
    a, c, e, f = coefficients(x, y, d)

    # Compute the integral
    α = d * sum(a)
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
function integrate(x::AbstractVector, y::AbstractVector, g::Function; d::Real=0.001)
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
