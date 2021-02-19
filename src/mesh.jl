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
function getpoint(ngon::Ngon; maxiter::Int=100_000) 
    tess = spt.Delaunay(coordinates(ngon))
    A, b = box(ngon) 
    iter = 1 
    while iter ≤ maxiter 
        val = A * rand(2) + b 
        # pnt = Point(val[1], val[2]) 
        pnt = Point(BigFloat(val[1]), BigFloat(val[2])) 
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

"""
    $SIGNATURES

Returns the two-dimensional projections of the three-dimensional points `pnts3d`
"""
project(pnts3d::AbstractVector{<:Point}) = [Point(pnt[1], pnt[2]) for pnt in pnts3d]

"""
    $SIGNATURES

Returns a two dimensional mesh whose points are the projections of `msh3`. 
"""
project(msh3::GeometryBasics.Mesh) = GeometryBasics.Mesh(project(msh3.position), copy(faces(msh3)))

"""
    $SIGNATURES

Returns a tuple of a Delaunay tessellation and a correspoding three dimensional mesh. 
"""
function triangulate(pnts3d) 
    pnts2d = project(pnts3d)
    tess = spt.Delaunay(pnts2d)
    trifaces = [TriangleFace(val...) for val in eachrow(tess.simplices .+ 1)]
    tess, GeometryBasics.Mesh(pnts3d, trifaces)
end 

"""
    $SIGNATURES 

Returns the bounding triangle of the points `pnts3d`.
"""
function findouttriangle(pnts3d)
    pnts2d = project(pnts3d)
    # TODO: Call c code directly
    hull = spt.ConvexHull(collect(hcat(collect.(pnts2d)...)') )
    if length(hull.vertices) == 3
        # If the convex hull is a triangle, just return the triangle 
        Triangle(Point.(pnts3d[hull.vertices .+ 1])...)
    else 
        # If the the convex hull is not a triangle but a polygon, construct the 
        # boundary polygon, construct a mesh from the polygon and return the 
        # maximum triangle with the maximum area.
        polygon = Ngon(
            SVector{length(hull.vertices)}(Point.(pnts3d[hull.vertices .+ 1]))
        )
        msh = GeometryBasics.mesh(polygon)
        idx = argmax([area(trig.points) for trig in msh])
        msh[idx]
    end 
end

"""
    $SIGNATURES

Returns a fractal surface interpolation function that interpolates `pnts3d`.
"""
function interpolate(pnts3d; α=0.01, f0 = (x, y) -> 0., maxiters=15, gettransforms::Bool=false)
    # Fint convex hull, i.e, the boundary triangle 
    outtrig = findouttriangle(pnts3d)

    # Consruct a 3d mesh 
    tess, msh3 = triangulate(pnts3d)

    # Compute transforms 
    transforms = map(intrig -> gettransform(outtrig, intrig, α), msh3)

    # Define mapping 
    # footprints = Point2{Float64}[]      # Uncomment to make a container for the footprints of any point 
    function mapping(f)
		function fnext(xd, yd)
            pnt = Point(xd, yd)
            # push!(footprints, pnt)
			idx = tess.find_simplex(pnt)[1] + 1  
            idx == 0 && return NaN 
            transform = transforms[idx]
            A, b = transform.A, transform.b
			linv = A[1:2, 1:2] \ (pnt - b[1:2])
			A[3, 1 : 2] ⋅ linv + α * f(linv[1], linv[2]) + b[3]
        end
	end 

    # Main iteration of the fractal interpolation.
    interpolant = ∘((mapping for i in 1 : maxiters)...)(f0)
	gettransforms ? (interpolant, transforms) : interpolant
	# gettransforms ? (interpolant, transforms, footprints) : (interpolant, footprints)
end 

"""
    $SIGNATURES

Returns a named tuple of `A` and` `b` such that `L(x) = A * x + b` maps `outtrig` to `intrig`. 
"""
function gettransform(outtrig, intrig, α::Real=1.) 
    outmat = collect(hcat(coordinates(outtrig)...)')
    inmat = collect(hcat(coordinates(intrig)...)')
    inmat[:, end] -= α * outmat[:, end]
    outmat[:, end] = ones(3)
    sol = outmat \ inmat
    (A = [  sol[1, 1]   sol[2, 1]   0;
            sol[1, 2]   sol[2, 2]   0; 
            sol[1, 3]   sol[2, 3]   α           ],
    b = [   sol[3, 1],  sol[3, 2],  sol[3, 3]   ])
end 


# ---------------------------------- Plot recipe ------------------------------------------------------ # 


@recipe(Trisurf, msh, f) do scene 
    Attributes(
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
    mesh!(plt, msh3, color=plt.meshcolor3, colormap=plt.colormap, visible=plt.visible) 
    wireframe!(plt, msh3, linewidth=plt.wflinewidth3)
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

