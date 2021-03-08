# This file includes an example file for 1D interpolation. 

using FractalTools 
using Plots 

# Generate data 
xi, dx, xf = 0., 0.1, 10.
x = collect(xi : dx : xf) 
y = sin.(x) 
pts = collect.(zip(x, y))
interp = interpolate(pts, Interp1D(0.01));

# Calculate interpolant 
xt = collect(xi : 0.1dx : xf)
yt = interp.(xt)

# Plot interpolation 
plt = plot(xt, yt)
scatter!(plt, x, y, markersize=3)
