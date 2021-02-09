
using FractalTools 
using Makie 

# Construct tessellation 
p1, p2, p3 = [0., 0.], [1., 0.], [0.5, 1.]
dlntess = DelaunayTessellation(p1, p2, p3, addboundarypoints=true)

tess = Node(dlntess.tessellation)
apnt = Node([0., 0])
tpnts = Node(p1 .+ eps())
msh = @lift(tomesh($tess))

fig = Figure(resolution = (500, 500))
ax = fig[1, 1] = Axis(fig, resolution = (500, 500))
btn = Button(fig[2,1], tellwidth = false) 
hidedecorations!(ax)
mesh!(msh, color=:skyblue)
wireframe!(msh, color=:black, linewidth=3)

scene = fig.scene

on(scene.events.mouseposition) do val 
    pos = to_world(scene, Point2f0(val))
    @show pos
end 

on(scene.events.mousebuttons) do buttons
   if ispressed(scene, Mouse.left)
        pos = to_world(scene, Point2(scene.events.mouseposition[]))
        @info "$pos added"
        addpoint!(dlntess, pos)
        tess[] = dlntess.tessellation
        apnt[] = pos
        update!(scene)
   end
   return
end

display(fig)

# # Test point 
# tpnt = (p1 + p2) / 2 .+ 1
# tidx = simplices(dlntess)[locate(dlntess, tpnt), :]
# foo = points(dlntess)[tidx, :] 
# scatter!(scene, foo[:, 1], foo[:, 2], color=:red, markersize=10)
# scatter!(scene, [tpnt[1]], [tpnt[2]], color=:orange, markersize=10)

