
module Foo 

# TODO: Add refine(tess) function to remove excessice points 

using PyCall
using GeometryBasics
using Clustering 

export trimesh, clustercenters

spt = pyimport_conda("scipy.spatial", "scipy") 

function box(p1::AbstractVector, p2::AbstractVector, p3::AbstractVector)
    x, y = first.([p1, p2, p3]), last.([p1, p2, p3]) 
    xmin, xmax, ymin, ymax = minimum(x), maximum(x), minimum(y), maximum(y)
    xwidth, ywidth = xmax - xmin, ymax - ymin 
    A = [xwidth 0; 0 ywidth]
    b = [xmin, ymin] 
    A, b
end

struct DelaunayTessellation{TS,TR}
    tessellation::TS
    triangle::TR
end 
function DelaunayTessellation(p1::AbstractVector, p2::AbstractVector, p3::AbstractVector) 
    p1, p2, p3, p4 = Point(p1...), Point(p2...), Point(p3...), Point(((p1 + p2 + p3) / 3)...)
    triangle = Triangle(p1, p2, p3)
    pnts = append!([p1, p2, p3, p4],  boundarypoints(p1, p2, p3))
    tessellation = spt.Delaunay(pnts, incremental=true)
    DelaunayTessellation(tessellation, triangle)
end

function addpoint!(tess::DelaunayTessellation)
    tri = tess.triangle 
    A, b = box(tri.points...)
    val = A * rand(2) + b 
    pnt = Point(val[1], val[2]) 
    pnt ∈ tri && tess.tessellation.add_points([pnt])
    pnt 
end 

function clustercenters(pnts, numpoints) 
    mat = [getindex.(pnts, 1) getindex.(pnts, 2)]
    [Point(val[1], val[2]) for val in eachcol(kmeans(mat', numpoints).centers)]
end 

function boundarypoints(p1, p2, p3; numpoints::Int=10)
    vcat([p1, p2, p3], 
    linepoint(p1, p2, numpoints=numpoints), 
    linepoint(p1, p3, numpoints=numpoints), 
    linepoint(p2, p3, numpoints=numpoints))
end 

function linepoint(p1, p2; numpoints::Int=10) 
    x = collect(range(p1[1], p2[1], length=numpoints))
    y = collect(range(p1[2], p2[2], length=numpoints))
    [Point(xi, yi) for (xi, yi) in zip(x, y)]
end

function tomesh(tri)
    vs = [Point(val[1], val[2]) for val in eachrow(tri.points)]
    fs = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tri.simplices .+ 1)]
    Mesh(vs, fs)
end 

end # module 

import .Foo 
p1, p2, p3 = [0., 0.], [1., 0.], [0.5, 1.]
tess = Foo.DelaunayTessellation(p1, p2, p3) 
pnt = Foo.addpoint!(tess)
msh = Foo.tomesh(tess.tessellation)
import Makie 
fig, ax, plt = Makie.mesh(msh, color=last.(msh.position))
Makie.wireframe!(msh) 
Makie.scatter!(getindex.(msh.position, 1), getindex.(msh.position, 2))
Makie.scatter!(getindex.([p1, p2, p3], 1), getindex.([p1, p2, p3], 2), color=:red)
display(fig) 
