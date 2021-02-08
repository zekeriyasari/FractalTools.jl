using VoronoiDelaunay 
using Makie 
using GeometryBasics 

const wdt = max_coord - min_coord 

function addpoint!(tess) 
    p = VoronoiDelaunay.Point(min_coord + rand() * wdt, min_coord + rand() * wdt)
    push!(tess, p) 
    tess
end

function plotedge(ax, edge) 
    x = [getx(geta(edge)), getx(getb(edge))]
    y = [gety(geta(edge)), gety(getb(edge))]
    lines!(ax, x, y, linewidth=3) 
    scatter!(ax, x, y, markersize=10, color=:black) 
end 

function getmesh(trigs)
    i = 1 
    pnts = GeometryBasics.Point2{Float64}[]
    fcs = TriangleFace{Int64}[]
    for trig in trigs 
        push!(pnts, GeometryBasics.Point(getx(trig._a), gety(trig._a)))
        push!(pnts, GeometryBasics.Point(getx(trig._b), gety(trig._b)))
        push!(pnts, GeometryBasics.Point(getx(trig._c), gety(trig._c)))
        push!(fcs, TriangleFace(i, i + 1, i + 2))
        i += 3
    end 
    pnts, fcs
end 

tess = DelaunayTessellation()
foreach(i -> addpoint!(tess), 1 : 10) 
trigs = filter(!isexternal, tess._trigs)
ps, fs = getmesh(trigs)
