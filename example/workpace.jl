# Test file 

using FractalTools
using Plots
pyplot()
default(:size, (1900, 1000))

# Construct interpolation data 
xi, dx, xf = 0, 0.1, 1
x = collect(xi : dx : xf) 
y = sin.(2π * x)
d = 0.1 * ones(length(x) - 1)

plt = plot(layout=(5,1))
for niter in 1 : 5
    interp = fif(x, y, d, niter=niter)
    xd = collect(xi : 0.1dx : xf)
    plot!(xd, interp.(xd), label="$niter", subplot=niter)
    scatter!(x, y, subplot=niter)
end
display(plt)

