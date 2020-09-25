# Test file 

using FractalTools
using Plots 

# Construct interpolation data 
xi, dx, xf = 0, 0.1, 1
x = collect(xi : dx : xf) 
y = sin.(2π * x)
d = 0.01 * ones(length(x) - 1)

plt = plot()
for niter in 1 : 1
    interp = fif(x, y, d, niter=niter)
    xd = collect(xi : 0.1dx : xf)
    plot!(xd, interp.(xd))
end
display(plt)

