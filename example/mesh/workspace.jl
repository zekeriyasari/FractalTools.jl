using GeometryBasics 
using Makie 
using Random 
using FractalTools
import AbstractPlotting
Random.seed!(0)



# Define function 
pts = [
    Point(1.25, 1.25), 
    Point(1.75, 1.75), 
    Point(1.5, 1.75)
]
trig = Triangle(pts...)
f(x, y) = (x / 5)^2 + (y / 5)^2
pnts3d = getdata(f, trig, 100)
pnts2d = project(pnts3d)

tess, msh3 = triangulate(pnts3d)
msh2 = project(msh3) 

fig, ax, plt = trisurf(msh2, f)
wireframe!(msh2)

itp = interpolate(pnts3d)


pnt = pts[2]
itp(pnt[1], pnt[2]) 

ppnt = Point(-1.477162590008043e-14, 2.0000000000000107)
scatter!(pnt, color=:green, markersize=20)
scatter!(ppnt, color=:red, markersize=20)

display(fig) 

