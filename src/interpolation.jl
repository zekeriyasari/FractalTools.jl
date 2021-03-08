# This includes interpolation methods 

const PointVector{Dim} = AbstractVector{<:AbstractPoint{Dim, T}} 

abstract type CurveInterp end
abstract type SurfaceInterp end

struct Interp1D{T<:Union{<:Real, <:AbstractVector{<:Real}}} <: CurveInterp
    freevars::T 
end 

struct Interp2D{T<:Union{<:Real, <:AbstractVector{<:Real}}} <: SurfaceInterp
    freevars::T
end 

struct HInterp1D{T<:Union{<:AbstractMatrix, <:AbstractVector{<:AbstractMatrix}}} <: CurveInterp
    freevars::T 
end 

struct HInterp2D{T<:Union{<:AbstractMatrix, <:AbstractVector{<:AbstractMatrix}}} <: SurfaceInterp
    freevars::T 
end 


struct Tessellation{T<:Union{<:LineString, <:PyObject}}
    tess::T 
end 


function interpolate(pts::PointVector{2}, method::Interp1D; f0 = x -> 0., niter::Int = 10) 
    tess = tessellate(domainpts, method)
    domains = partition(pts, method) 
end 

function interpolate(pts::PointVector{3}, method::Interp2D; f0 = (x, y) -> 0., niter::Int = 10) 
    tess = tessellate(domainpts, method)
    domains = partition(pts, method) 
end 

function interpolate(pts::PointVector{3}, method::HInterp1D; f0 = x -> [0., 0.], niter::Int = 10) 
    tess = tessellate(domainpts, method)
    domains = partition(pts, method) 
end 

function interpolate(pts::PointVector{4}, method::HInterp2D; f0 = (x, y) -> [0., 0.], niter::Int = 10) 
    tess = tessellate(domainpts, method)
    domains = partition(pts, method) 
end 


project(pts::AbstractPoint, drop::Int=1) = [Point(pnt[1 : end - drop]...) for pnt in pts]
project(pts::AbstractPoint{2}, ::Interp1D)   = project(pts, 1)
project(pts::AbstractPoint{3}, ::HInterp1D)  = project(pts, 2)
project(pts::AbstractPoint{3}, ::Interp2D)   = project(pts, 1)
project(pts::AbstractPoint{4}, ::HInterp2D)  = project(pts, 2)


tessellate(pts::PointVector{2}, ::Interp1D)  = Tessellation(LineString(project(pts, method)))     
tessellate(pts::PointVector{3}, ::HInterp1D) = Tessellation(LineString(project(pts, method)))     
tessellate(pts::PointVector{3}, ::Interp2D)  = Tessellation(spt.Delaunay(project(pts, method)))
tessellate(pts::PointVector{4}, ::HInterp2D) = Tessellation(spt.Delaunay(project(pts, method)))


partition(pts::PointVector, method::CurveInterp) = LineString(pts) 
function partition(pts::PointVector, method::SurfaceInterp, domaintess::Tessellation=tessellate(pts, method)) 
    trifaces = [TriangleFace(val...) for val in eachrow(domaintess.tess.simplices .+ 1)]
    GeometryBasics.Mesh(pts, trifaces)
end 


getboundary(pts::PointVector, ::CurveInterp) = Line(pts[1], pts[end])
function getboundary(pts::PointVector, ::SurfaceInterp) 
    hull = spt.ConvexHull(collect(hcat(collect.(project(pts))...)') )
    if length(hull.vertices) == 3
        # If the convex hull is a triangle, just return the triangle 
        Triangle(Point.(pts[hull.vertices .+ 1])...)
    else 
        # If the the convex hull is not a triangle but a polygon, construct the 
        # boundary polygon, construct a mesh from the polygon and return the 
        # maximum triangle with the maximum area.
        polygon = Ngon(
            SVector{length(hull.vertices)}(Point.(pts[hull.vertices .+ 1]))
        )
        msh = GeometryBasics.mesh(polygon)
        idx = argmax([area(trig.points) for trig in msh])
        msh[idx]
    end 
end 


# Interp1D
function gettransform(outline::Line, inline::Line, freevar::Real) 
    outmat = collect(hcat(coordinates(outline)...)')
    inmat = collect(hcat(coordinates(inline)...)')
    inmat[:, end] -= outmat[:, end] * freevar
    outmat[:, end] .= 1
    a11, b1, a21, b2 = outmat \ inmat
    A = [a11    0; 
         a21    freevar]
    b = [b1, b2]
	Transformation(A, b)
end 

# HInterp1D
function gettransform(outline::Line, inline::Line, freevar::AbstractMatrix) 
    outmat = collect(hcat(coordinates(outline)...)')
    inmat = collect(hcat(coordinates(inline)...)')
    inmat[:, 2 : 3] -= outmat[:, 2 : 3] * freevar'
    outmat = [outmat[:, 1] ones(2)]
    a11, b1, a21, b2, a31, b3 = outmat \ inmat
    A = [a11    0               0; 
         a21    freevar[1,1]    freevar[1,2];
         a31    freevar[2,1]    freevar[2,2]
         ]
    b = [b1, b2, b3]
	Transformation(A, b)
end 

# Interp2D
function gettransform(outtrig::Triangle, intrig::Triangle, freevar::Real) 
    outmat = collect(hcat(coordinates(outtrig)...)')
    inmat = collect(hcat(coordinates(intrig)...)')
    inmat[:, end] -= outmat[:, end] * freevar
    outmat[:, end] .= 1
    a11, a12, b1, a21, a22, b2, a31, a32, b3 = outmat \ inmat
    A = [a11    a12    0; 
         a21    a22    0;
         a31    a32    freevar
         ]
    b = [b1, b2, b3]
	Transformation(A, b)
end 

# HInterp2D
function gettransform(outtrig::Triangle, intrig::Triangle, freevar::AbstractMatrix) 
    outmat = collect(hcat(coordinates(outtrig)...)')
    inmat = collect(hcat(coordinates(intrig)...)')
    inmat[:, 3:4] -= outmat[:, 3:4] .* freevar'
    outmat = [outmat[:, 1:2] ones(3)]
    a11, a12, b1, a21, a22, b2, a31, a32, b3, a41, a42, b4 = outmat \ inmat
    A = [a11    a12    0                0; 
         a21    a22    0                0;
         a31    a32    freevar[1, 1]    freevar[1, 2];
         a41    a42    freevar[2, 1]    freevar[2, 2];
         ]
    b = [b1, b2, b3, b4]
	Transformation(A, b)
end 
