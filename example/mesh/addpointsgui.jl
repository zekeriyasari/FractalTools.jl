using FractalTools 
using Makie 

# Contruct a tesslation 
dlntess = Node(DelaunayTessellation([0., 0.], [1., 0], [0.5, 1.], addboundarypoints=false))

# Definea an observables
msh = lift(dlntess) do val 
    @info "change in dlntess"
    tomesh(val.tessellation)
end 
vtx = lift(msh) do val 
    @info "change in msh"
    val.position 
end 

# Plot tesselation 
scene = Scene(resolution=(500, 500), scale_plot = false)
mesh!(scene, msh, color=:lightblue)
wireframe!(scene, msh, color=:black)
scatter!(scene, vtx, color=:black)

# Add mouseposition interaction 
on(scene.events.mouseposition) do pos 
    val = to_world(scene, Point(pos))
    @show val 
end 

# Add mouse button interaction 
on(scene.events.mousebuttons) do buttons
   if ispressed(scene, Mouse.left)
        pos = to_world(scene, Point(scene.events.mouseposition[]))
        @info "$pos added"
        addpoint!(dlntess[], pos)
        dlntess[] = dlntess[]   # Trigger change
        push!(vtx[], pos)
   end
end


# Display scene 
display(scene)
