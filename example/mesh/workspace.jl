
using FractalTools 
using Makie 
using AbstractPlotting
import AbstractPlotting.GeometryBasics: TriangleFace 

# Construct tessellation 
p1, p2, p3 = [-1., √3.], [-1., -√3], [2., 0.]
dlntess = DelaunayTessellation(p1, p2, p3)

while npoints(dlntess) ≤ 1000 
    addpoint!(dlntess) 
end 

finegrain!(dlntess, 100)

# Define function 
f(x, y) = x^2 + y^2

# Plot the result 
fig, ax, plt = tessplotf(dlntess, f, vmarkersize=10)
wireframe!(tomesh(dlntess.tessellation), linewidth=3)
display(fig)
