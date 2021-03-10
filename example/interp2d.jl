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

function err(x, y)
    fval = f(x, y) 
    ival = interp(x, y) 
    abs(fval - ival) / abs(fval) * 100
end

tpts = getdata(ngon, npts)

fvals = map(pnt -> Point(f(pnt...)...), tpts)
ivals = map(pnt -> Point(interp(pnt...)...), tpts)


fig, ax, plt = trisurf(tpts, f, meshcolor3=:red)
trisurf!(ax, tpts, interp, meshcolor3=:blue)
trisurf(tpts, evals, meshcolor3=last.(evals), colormap=:viridis)

