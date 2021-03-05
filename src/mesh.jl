#  This file includes mesh related methods. 

export boundarypoints, isvalidpoint, getpoint, box, project, triangulate, findouttriangle, interpolate, gettransform, disperse, getdata 
export trisurf 

"""
    $SIGNATURES 

Returns boundary points of on the edges of the triangle `p1, p2, p3`. `numpoints` is the number of points on the edges of the triangle.
"""
function boundarypoints(ngon::Ngon; numpoints::Int=10)
    vcat(
        map(((pnt1, pnt2),) -> linepoint(pnt1, pnt2), TupleView{2, 1}(SVector([ngon.points; [ngon.points[1]]]...)))...
        )
    # vcat(
    #     linepoint(p1, p2, numpoints=numpoints), 
    #     linepoint(p1, p3, numpoints=numpoints), 
    #     linepoint(p2, p3, numpoints=numpoints)
    #     )
end 

"""
    $SIGNATURES

Returns the points of the edge whose end points are `p1` and `p2`. `numpoints` is the number of points. 
"""
function linepoint(p1::Point, p2::Point; numpoints::Int=10) 
    x = collect(LinRange(p1[1], p2[1], numpoints))
    y = collect(LinRange(p1[2], p2[2], numpoints))
    [Point(xi, yi) for (xi, yi) in zip(x, y)]
end

""""
    $SIGNATURES

Returns true if `pnt` is in `tess`. 
"""
isvalidpoint(pnt::Point, tess) = tess.find_simplex(pnt)[1] ≥ 0 && pnt !== Point(NaN, NaN)


"""
    $SIGNATURES

Returns a random valid point in `ngon`. `maxiters` is the number of iteration while finding the point. 
"""
function getpoint(ngon::AbstractPolygon{Dim, T}; maxiter::Int=100_000) where {Dim, T}
    tess = spt.Delaunay(coordinates(ngon))
    A, b = box(ngon) 
    iter = 1 
    while iter ≤ maxiter 
        val = A * rand(T, 2) + b 
        pnt = Point(val[1], val[2]) 
        isvalidpoint(pnt, tess) && return pnt
        iter += 1
    end 
    return Point(NaN, NaN)  # For type-stability 
end 

"""
    $SIGNATURES 

Returns a vector of random points that are dispersed into `ngon`. `npoints` is the number of points to be dispersed. 
"""
function disperse(ngon, npoints) 
    allpnts = [getpoint(ngon) for i in 1 : 10 * npoints]
    ctrpnts = [Point(val[1], val[2]) for val in eachcol(kmeans(hcat(collect.(allpnts)...), npoints).centers)]
    vcat(ngon.points, boundarypoints(ngon), ctrpnts)
end

"""
    $SIGNATURES

Filters under quality triangles
"""
function filtertriangle(quality, args...; kwargs...) end 

"""
    $SIGNATURES

Returns a three-dimensional interpolation data `pnts`. `pnts` is a vector of three-dimensional points `pi = Point(xi, yi, zi)` where `xi` and `yi` are from the dispersed points and `zi = f(xi, yi)`. 
"""
getdata(f, ngon::Ngon, npts::Int) = [Point(pnt[1], pnt[2], f(pnt[1], pnt[2])) for pnt in disperse(ngon, npts)] 


"""
    $SIGNATURES 

Returns the matrix `A` and a vector `b` such that for any point `p` the transformation `T(p) = A * p + b` moves `p` into the bounding box of the ngon formed by the points `p1, p2, p3`.  
"""
function box(ngon::Ngon)
    coords = coordinates(ngon) 
    x, y = getindex.(coords, 1), getindex.(coords, 2) 
    xmin, xmax, ymin, ymax = minimum(x), maximum(x), minimum(y), maximum(y)
    xwidth, ywidth = xmax - xmin, ymax - ymin 
    A = [xwidth 0; 0 ywidth]
    b = [xmin, ymin] 
    A, b
end


# ---------------------------------- Plot recipe ------------------------------------------------------ # 


@recipe(Trisurf, msh, f) do scene 
    AbstractPlotting.Attributes(
        wireframe2 = false,
        wfcolor = :black,
        wflinewidth = 2,
        vmarkercolor = :red,
        wflinewidth3 = 3,
        vmarkercolor3 = :orange, 
        vmarkersize3 = 20,
        meshcolor3 = :black,
        colormap = :viridis,
        visible = true
    )
end

function AbstractPlotting.plot!(plt::Trisurf) 
    msh3 = plt[1][]
    AbstractPlotting.mesh!(plt, msh3, color=plt.meshcolor3, colormap=plt.colormap, visible=plt.visible) 
    AbstractPlotting.wireframe!(plt, msh3, linewidth=plt.wflinewidth3)
    plt
end

function AbstractPlotting.convert_arguments(::Type{<:Trisurf}, msh2::GeometryBasics.Mesh, f::Function)
    msh3 = GeometryBasics.Mesh([Point(pnt[1], pnt[2], f(pnt[1], pnt[2])) for pnt in msh2.position], faces(msh2))
    return (msh3, msh2, f)
end 

function AbstractPlotting.convert_arguments(::Type{<:Trisurf}, pnts::AbstractVector{<:Point2}, f::Function) 
    tess = spt.Delaunay(pnts) 
    _position = [Point(pnt[1], pnt[2], f(pnt[1], pnt[2])) for pnt in eachrow(tess.points)]
    _faces = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.simplices .+ 1)]
    msh3 = GeometryBasics.Mesh(_position, _faces)
    return (msh3, pnts, f)
end 

function AbstractPlotting.convert_arguments( ::Type{<:Trisurf}, pnts2d::AbstractVector{<:Point2}, pnts1d::AbstractVector{<:Point1})
    tess = spt.Delaunay(pnts2d) 
    _position = [Point(pnt2[1], pnt2[2], pnt1[1]) for (pnt2, pnt1) in zip(pnts2d, pnts1d)]
    _faces = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.simplices .+ 1)]
    msh3 = GeometryBasics.Mesh(_position, _faces)
    return (msh3, pnts2d, pnts1d)
end 

function AbstractPlotting.convert_arguments( ::Type{<:Trisurf}, pnts3d::AbstractVector{<:Point3})
    pnts2d = project(pnts3d)
    tess = spt.Delaunay(pnts2d) 
    fcs = [TriangleFace(val[1], val[2], val[3]) for val in eachrow(tess.simplices .+ 1)]
    msh3 = GeometryBasics.Mesh(pnts3d, fcs)
    return (msh3, pnts3d)
end 

