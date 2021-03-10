# This file includes and example file for 2D interpolation 

using FractalTools 
using GeometryBasics
using Makie 

# Generate data 
f(x, y) = x^2 + y^2 + 1
ngon = Triangle(
    Point(BigFloat(0.), BigFloat(0.)), 
    Point(BigFloat(1.), BigFloat(0.)), 
    Point(BigFloat(0.5), BigFloat(1.)))
npts = 100
pts = getdata(f, ngon, npts)
interp = interpolate(pts, Interp2D(0.01))

tpts = getdata(ngon, npts)

fvals = map(pnt -> Point(f(pnt...)...), tpts)
ivals = map(pnt -> Point(interp(pnt...)...), tpts)
evals = (fvals - ivals) ./ abs.(fvals) * 100)
@show maximum(evals)


fig, ax, plt = trisurf(tpts, ivals, meshcolor3=last.(fvals), colormap=:viridis)
trisurf!(ax, tpts, ivals, meshcolor3=last.(ivals), colormap=:heat)
trisurf(tpts, evals, meshcolor3=last.(evals), colormap=:viridis)

