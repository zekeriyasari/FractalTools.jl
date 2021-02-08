
module Foo 

using VoronoiDelaunay 
using GeometricalPredicates

export Domain, addpoint!

const minc, maxc = VoronoiDelaunay.min_coord, VoronoiDelaunay.max_coord

struct Domain{T, S}
    boundary::T 
    tessellation::S 
    function Domain{T, S}(boundary, tessellation) where {T, S} 
        foreach(checkpoint, [boundary._a, boundary._b, boundary._c])
        foreach(p -> push!(tessellation, p), [boundary._a, boundary._b, boundary._c])
        new{T, S}(boundary, tessellation)
    end
end 
Domain(boundary::T, tesselation::S) where {T, S} = Domain{T, S}(boundary, tesselation)
Domain(p1::Point2D, p2::Point2D, p3::Point2D) = (foreach(checkpoint, [p1, p2, p3]); Domain(Triangle(p1, p2, p3), DelaunayTessellation()))
Domain(pts::Real...) = Domain(Triangle(pts...), DelaunayTessellation())
Domain() = Domain(minc, minc, maxc, minc, (maxc + minc) / 2, maxc)

function checkpoint(p::Point2D) 
    ((minc ≤ p._x ≤ maxc) && (minc ≤ p._y ≤ maxc)) || 
        error("Expected coordinate in the inteval [$minc, $maxc], got $p") 
end


function addpoint!(domain::Domain)
    val = rand(2) + [minc, minc]
    p = Point(val[1], val[2]) 
    intriangle(domain.boundary, p) == 1 && push!(domain.tessellation, p)
    domain
end

end # module 

using .Foo 
using VoronoiDelaunay 
using Makie 

function filldomain(dm)
    idx = 0
    while dm.tessellation._total_points_added ≤ 10
        Foo.addpoint!(dm)
        idx += 1
        idx ≥ 1000 && break 
    end 
end 

dm = Node(Domain())
x = @lift(getplotxy(delaunayedges($dm.tessellation))[1])
y = @lift(getplotxy(delaunayedges($dm.tessellation))[2])
fig = Figure() 
ax = Axis(fig[1, 1])
lines!(ax, x, y)
scatter!(ax, x, y)

btn = Button(fig[2, 1], tellwidth=false, label="addpoint")
on(btn.clicks) do n
    @show n  
    addpoint!(dm[])
    @show filter(!isexternal, dm[].tessellation._trigs)
end

