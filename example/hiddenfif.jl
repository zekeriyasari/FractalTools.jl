# This file includes one dimensional interpolation. 

using FractalTools 
using Plots 

# Generate interpolation data 
xi, dx, xf = 0, 0.1, 1
x = collect(xi : dx : xf) 
y = sin.(2π * x)
z = cos.(2π * x)
d = 0.1 * ones(length(x) - 1)
h = 0.1 * ones(length(x) - 1)
l = 0.1 * ones(length(x) - 1)
m = 0.1 * ones(length(x) - 1)

# Construct interpolant 
# ifs = gethiddenifs(x, y, z, d, h, l, m)
interp = hiddenfif(x, y, z, d, h, l, m)

# Plot data 
xd = collect(xi : 0.01dx : xf)
yd = interp.(xd)
plot(xd, getindex.(yd, 1))
scatter!(x, y, marker=(:circle, 2))

