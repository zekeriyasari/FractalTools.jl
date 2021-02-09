# This file 

export DelaunayTessellation, addpoint!, tomesh, npoints, locate, finegrain!, points, simplices, getpoint, tessplot, tessplotf

"""
    $TYPEDEF

Delaunay tessellation 
"""
mutable struct DelaunayTessellation{TS,TR}
    "Tessellation"
    tessellation::TS
    "Boundary triangle"
    triangle::TR
end 

"""
    $SIGNATURES 

Construct a Delaunay triangulation from the points `p1, p2, p3`
"""
function DelaunayTessellation(p1::AbstractVector, p2::AbstractVector, p3::AbstractVector; addboundarypoints::Bool=true) 
    p1, p2, p3, p4 = Point(p1...), Point(p2...), Point(p3...), Point(((p1 + p2 + p3) / 3)...)
    triangle = Triangle(p1, p2, p3)
    pnts = [p1, p2, p3, p4]
    addboundarypoints && append!(pnts,  boundarypoints(p1, p2, p3))
    tessellation = spt.Delaunay(pnts, incremental=true)
    DelaunayTessellation(tessellation, triangle)
end

"""
    $SIGNATURES

Adds a point inside the tessellation by preserving delauneyhood.
"""
function addpoint!(dlntess::DelaunayTessellation, pnt::Point2=getpoint(dlntess))
    isvalidpoint(pnt, dlntess) && dlntess.tessellation.add_points([pnt])
    pnt 
end 

isvalidpoint(pnt::AbstractVector, dlntess::DelaunayTessellation) = pnt ∈ dlntess.triangle && pnt !== Point(NaN, NaN)

# Returns a valid point inside the triangle of the dlntess
function getpoint(dlntess::DelaunayTessellation; maxiter::Int=100_000) 
    tri = dlntess.triangle 
    A, b = box(tri.points...) 
    iter = 1 
    while iter ≤ maxiter 
        val = A * rand(2) + b 
        pnt = Point(val[1], val[2]) 
        isvalidpoint(pnt, dlntess) && return pnt
        iter += 1
    end 
    return Point(NaN, NaN)  # For type-stability 
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


"""
    $SIGNATURES 

Constructs a GeometryBasic.Mesh from `dlntess` ready for plotting.
"""
function tomesh(tri)
    vs = [Point(val[1], val[2]) for val in eachrow(tri.points)]
    fs = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tri.simplices .+ 1)]
    GeometryBasics.Mesh(vs, fs)
end 

# Returns the bouding box to the tessellation 
function box(p1::AbstractVector, p2::AbstractVector, p3::AbstractVector)
    x, y = first.([p1, p2, p3]), last.([p1, p2, p3]) 
    xmin, xmax, ymin, ymax = minimum(x), maximum(x), minimum(y), maximum(y)
    xwidth, ywidth = xmax - xmin, ymax - ymin 
    A = [xwidth 0; 0 ywidth]
    b = [xmin, ymin] 
    A, b
end


"""
    $SIGNATURES 

Returns the number of points in `dlntess`.
"""
npoints(dlntess::DelaunayTessellation) = dlntess.tessellation.npoints

"""
    $SIGNATURES 

Returns the coordinates of the points of `dlntess`. 
"""
points(dlntess::DelaunayTessellation) = dlntess.tessellation.points

"""
    $SIGNATURES 

Returns the simplices of `dlntess` 
"""
simplices(dlntess::DelaunayTessellation) = dlntess.tessellation.simplices .+ 1

"""
    $SIGNATURES 

Returns the triangle index of `dlntess` in which the `point` is.
"""
function locate(dlntess::DelaunayTessellation, point::AbstractVector)
    idx = dlntess.tessellation.find_simplex(point)[1] + 1
    idx == 0 ? error("The point $point cannot be located") : idx
end 

function finegrain!(dlntess::DelaunayTessellation, npnts::Int)   
    pnts = [Point(val[1], val[2]) for val in eachrow(dlntess.tessellation.points)]
    centers = clustercenters(pnts, npnts) 
    dlntess.tessellation = spt.Delaunay(append!(centers, boundarypoints(dlntess.triangle.points...)), incremental=true)
    dlntess
end 

function clustercenters(pnts, numpoints) 
    mat = [getindex.(pnts, 1) getindex.(pnts, 2)]
    [Point(val[1], val[2]) for val in eachcol(kmeans(mat', numpoints).centers)]
end 

####  Define a plotting recipe 

@recipe(TessPlot, dlntess) do scene 
    AbstractPlotting.Attributes(
        vcolor = :red,
        vmarkersize = 10,
        lwidth = 3
    )
end 

function AbstractPlotting.plot!(tessplot::TessPlot)
    msh = tomesh(tessplot[:dlntess][].tessellation)
    AbstractPlotting.mesh!(tessplot, msh, color=first.(msh.position))
    AbstractPlotting.wireframe!(tessplot, msh, linewidth=tessplot.lwidth) 
    AbstractPlotting.scatter!(tessplot, getindex.(msh.position, 1), getindex.(msh.position, 2), 
        color=tessplot.vcolor, markersize=tessplot.vmarkersize)
    tessplot
end 

"""
    tessplot(dlntess::DelaunayTessellation)

Plots the DelaunayTessellation `dlntess`
"""
function tessplot end 


@recipe(TessPlotF, dlntess, f) do scene 
    AbstractPlotting.Attributes(
        vcolor = :red,
        vmarkersize = 10,
        lwidth = 3
    )
end 

function AbstractPlotting.plot!(tessplotf::TessPlotF)
    _dlntess = tessplotf[:dlntess][] 
    _f = tessplotf[:f][]
    pnts = points(_dlntess)
    trigs = simplices(_dlntess)
    pnts3d = [Point3f0(pnt..., _f(pnt...)) for pnt in eachrow(pnts)]
    trigs3d = [TriangleFace(trig...) for trig in eachrow(trigs)]
    msh = GeometryBasics.Mesh(pnts3d, trigs3d)
    AbstractPlotting.mesh!(tessplotf, msh, color=first.(msh.position))
    AbstractPlotting.wireframe!(tessplotf, msh, linewidth=tessplotf.lwidth) 
    AbstractPlotting.scatter!(tessplotf, 
        getindex.(msh.position, 1), getindex.(msh.position, 2), getindex.(msh.position, 3),
        color=tessplotf.vcolor, markersize=tessplotf.vmarkersize)
    tessplotf
end 

"""
    tessplotf(dlntess::DelaunayTessellation, f)

Plots the DelaunayTessellation `dlntess`
"""
function tessplotf end 


# """
#     tessplot(dlntess::DelaunayTessellation)

# Plots the DelaunayTessellation `dlntess`

#     tessplot(ax::Axis, dlntess::DelaunayTessellation) 

# Plots `dlntess` on `ax`.
# """
# function tessplot end 

# @recipe(TessPlot, msh) do scene 
#     AbstractPlotting.Attributes(
#         vcolor = :red,
#         vmarkersize = 10,
#         lwidth = 3
#     )
# end 

# function AbstractPlotting.plot!(tessplot::TessPlot)
#     _msh = tessplot[:msh]
#     coords = _msh[].position
#     AbstractPlotting.mesh!(tessplot, _msh, color=first.(coords))
#     AbstractPlotting.wireframe!(tessplot, _msh, linewidth=tessplot.lwidth) 
#     AbstractPlotting.scatter!(tessplot, getindex.(coords, 1), getindex.(coords, 2), 
#         color=tessplot.vcolor, markersize=tessplot.vmarkersize)
#     tessplot
# end 
