# This file includes methods for data generation. 

export getdata, getpoint

"""
    $SIGNATURES

Returns a three-dimensional interpolation data `pnts`. `pnts` is a vector of three-dimensional points `pi = Point(xi, yi,
zi)` where `xi` and `yi` are from the dispersed points and `zi = f(xi, yi)`. 
"""
getdata(f, ngon::Ngon, npts::Int) = [Point(pnt[1], pnt[2], f(pnt[1], pnt[2])) for pnt in disperse(ngon, npts)] 

"""
    $SIGNATURES

Returns a random valid point in `ngon`. `maxiters` is the number of iteration while finding the point. 
"""
function getpoint(ngon::AbstractPolygon{Dim, T}; maxiter::Int=100_000) where {Dim, T}
    tess = spt.Delaunay(coordinates(ngon))
    A, b = boundboxtransforms(ngon) 
    iter = 1 
    while iter ≤ maxiter 
        val = A * rand(T, 2) + b 
        pnt = Point(val[1], val[2]) 
        isvalidpoint(pnt, tess) && return pnt
        iter += 1
    end 
    return Point(NaN, NaN)  # For type-stability 
end 

boundarypoints(ngon::Ngon; numpoints::Int=10) = vcat(
        map(
            ((pnt1, pnt2),) -> linepoint(pnt1, pnt2), TupleView{2, 1}(SVector([ngon.points; [ngon.points[1]]]...))
            )...
        )

function linepoint(p1::AbstractPoint, p2::AbstractPoint; numpoints::Int=10) 
    x = collect(LinRange(p1[1], p2[1], numpoints))
    y = collect(LinRange(p1[2], p2[2], numpoints))
    [Point(xi, yi) for (xi, yi) in zip(x, y)]
end

isvalidpoint(pnt::AbstractPoint, tess) = tess.find_simplex(pnt)[1] ≥ 0 && pnt !== Point(NaN, NaN)

function disperse(ngon::Ngon, npoints::Int) 
    allpnts = [getpoint(ngon) for i in 1 : 10 * npoints]
    ctrpnts = [Point(val[1], val[2]) for val in eachcol(kmeans(hcat(collect.(allpnts)...), npoints).centers)]
    vcat(ngon.points, boundarypoints(ngon), ctrpnts)
end

function boundboxtransforms(ngon::Ngon)
    coords = coordinates(ngon) 
    x, y = getindex.(coords, 1), getindex.(coords, 2) 
    xmin, xmax, ymin, ymax = minimum(x), maximum(x), minimum(y), maximum(y)
    xwidth, ywidth = xmax - xmin, ymax - ymin 
    A = [xwidth 0; 0 ywidth]
    b = [xmin, ymin] 
    A, b
end

# TODO: Complete function
function filtertriangle(quality, args...; kwargs...) end 
