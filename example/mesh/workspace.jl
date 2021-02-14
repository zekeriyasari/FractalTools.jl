using Makie 
using GeometryBasics 
using PyCall 
using FractalTools
using StaticArrays 

spt = pyimport_conda("scipy.spatial", "scipy") 

# Construct interpolation data 
hexagon = GeometryBasics.Ngon(@SVector [Point(BigFloat(cos(θ)), BigFloat(sin(θ))) for θ in 0 : π / 3 : 2π - π / 3])
outertess = spt.Delaunay(coordinates(hexagon))
tfcs = [TriangleFace(val...) for val in eachrow(outertess.simplices .+ 1)]
outermsh2 = GeometryBasics.Mesh(coordinates(hexagon), tfcs)
f(x, y) = x^2 + y^2
pnts3d = vcat(map(trig -> getdata(f, trig,100), outermsh2)...)
pnts2d = project(pnts3d)

tpnt = pnts2d[1]

pnts3dd = map(trig -> pnts3d[findall(pnt -> pnt ∈ trig, pnts2d)], outermsh2)
itps = map(pnts -> interpolate(pnts), pnts3dd)

function itp(xd, yd)
    pnt = Point(xd, yd)
    idx = outertess.find_simplex(pnt)[1] + 1
    itps[idx](xd, yd)    
end

function erf(xd, yd)
    fval = f(xd, yd) 
    ival = itp(xd, yd) 
    abs(fval - ival) / abs(fval) * 100
end

ivals = map(pnt -> Point(pnt[1], pnt[2], itp(pnt...)), pnts2d)
fvals = map(pnt -> Point(pnt[1], pnt[2], f(pnt...)), pnts2d)
evals = map(pnt -> Point(pnt[1], pnt[2], erf(pnt...)), pnts2d)

fig = Figure() 
ax1 = LScene(fig[1, 1])
ax2 = LScene(fig[1, 2])
trisurf!(ax1, ivals, meshcolor3=:green)
trisurf!(ax1, fvals, meshcolor3=:red)
trisurf(evals, meshcolor3=:red)

display(fig)

# Plots 
fig, ax, plt = scatter(pnts3d)
scatter!(pnts3dd[1], color=:red)
scatter!(pnts3dd[2], color=:green)
wireframe!(outermsh2[1].points, color=:red)
wireframe!(outermsh2[2].points, color=:red)
scatter!(pnts2d)
wireframe!(outermsh2, linewidth=2)
scatter!(tpnt, color=:red, markersize=20)
mesh!(trig)
display(fig)