using FractalTools 
using Makie 

# Construct a tesselation 
pts = [0., 0.], [1., 0], [0.5, 1.] 
dlntess = DelaunayTessellation(pts...) 

# Plot the tessellation. 
tessplot(dlntess)
