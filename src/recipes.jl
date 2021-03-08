# This file includes plot recipes for Makie

export trisurf

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
