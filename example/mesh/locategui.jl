using FractalTools 
using Makie 

# Contruct a tesslation 
dlntess = Node(DelaunayTessellation([0., 0.], [1., 0], [0.5, 1.], addboundarypoints=true))

# Definea an observables
msh = lift(dlntess) do val 
    tomesh(val.tessellation)
end 
vtx = lift(msh) do val 
    val.position 
end 

# Test point 
tpnt = Node(getpoint(dlntess[]))
tripts = lift(tpnt) do val 
    tidx = simplices(dlntess[])[locate(dlntess[], val), :]
    map(val -> Point(val...), eachrow(points(dlntess[])[tidx, :]))
end 

# Plot tesselation 
fig = Figure(resolution=(500, 500), scale_plot = false, camera = AbstractPlotting.campixel!)
ax = fig[1, 1] =  Axis(fig) 
mesh!(ax, msh, color=:lightblue)
wireframe!(ax, msh, color=:black)
scatter!(ax, vtx, color=:black)
scatter!(ax, tripts, color=:red, markersize=10)
scatter!(ax, tpnt, color=:orange, markersize=10)

# Add mouse button interaction 
on(ax.scene.events.mousebuttons) do buttons
   if ispressed(ax.scene, Mouse.left)
        area = pixelarea(ax.scene)[]
        mp = ax.scene.events.mouseposition[]
        mp = Point2f0(mp) .- minimum(area)
        pos = to_world(ax.scene, mp)
        @info "$pos added"
        tpnt[] = pos
   end
end

# Add button interaction.
on(ax.scene.events.keyboardbuttons) do button 
    if ispressed(button, Keyboard.a)
        pnt = addpoint!(dlntess[])
        @info "$pnt is added"
    end 
    if ispressed(button, Keyboard.f)
        finegrain!(dlntess[], 10)
        @info "Tessellation is fine grained"
    end 
    dlntess[] = dlntess[]   # Trigger change
end 

# Display scene 
display(fig)
