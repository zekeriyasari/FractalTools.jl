using FractalTools
using Makie 
using GeometryBasics

# Construct a triangular mesh 
msh = trimesh([-1., √3], [-1., -√3], [2, 0.])


# Plot the mesh 
scn = Makie.mesh(msh, color=last.(msh.position))
Makie.wireframe!(msh)

coords = coordinates(msh)
trigs = faces(msh) 

idx = Node(1)
x = @lift(getindex.(coords[trigs[$idx]], 1) |> collect)
y = @lift(getindex.(coords[trigs[$idx]], 2) |> collect)
scatter!(x, y)
