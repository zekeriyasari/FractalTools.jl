using GeometryBasics 
using Makie 
using Random 
Random.seed!(0)

# Define function 
pts = [
    Point(-1., -1.), 
    Point(1., -1.), 
    Point(0., 2.)
] .+ Point(2, 3) 
trig = Triangle(pts...)
f(x, y) = x^2 + y^2 
fig, ax, plt = mesh(trig) 
