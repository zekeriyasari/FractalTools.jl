# This file includes one dimensional interpolation. 

using FractalTools 
using Plots 

# Generate interpolation data 
xi, dx, xf = 0, 0.1, 1
x = collect(xi : dx : xf) 
y = sin.(2π * x)
d = 0.0 * ones(length(x) - 1)

# Construct interpolant 
interp = fif(x, y, d)

# Plot data 
xd = collect(xi : 0.01dx : xf)
yd = interp.(xd)
plot(xd, yd)
scatter!(x, y, marker=(:circle, 5))

