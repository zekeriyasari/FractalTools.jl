# This file includes and example file for 2D interpolation 

using FractalTools 
using Plots 
using GeometryBasics

# Generate data 
f(x, y) = x^2 + y^2 + 1
ngon = Triangle(
    Point(BigFloat(0.), BigFloat(0.)), 
    Point(BigFloat(1.), BigFloat(0.)), 
    Point(BigFloat(0.5), BigFloat(1.)))
npts = 100
pts = getdata(f, ngon, npts)
interp = interpolate(pts, Interp2D(0.01));
