using FractalTools
using Makie 

# Construct a triangular mesh 
msh = trimesh([-1., √3], [-1., -√3], [2, 0.])

# Define function
f(x, y) = (x^2 + y^2) 

# Plot triangular surface plots
msh3 = spanmesh(msh, f)

# Plot the surface
scn = Makie.mesh(msh3, color=last.(msh3.position))
Makie.wireframe!(msh3) 
Makie.mesh!(msh, color=last.(msh.position))
Makie.wireframe!(msh)
display(scn)
