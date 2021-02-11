
using FractalTools 
using GeometryBasics
using Makie 

# Construct delaunay 
p1, p2, p3 = [-1., -1], [1., -1], [0., 1.]
tridln = TriDelaunay(p1, p2, p3)

while npoints(tridln) ≤ 1000 
    addpoint!(tridln) 
end 

finegrain!(tridln, 100)

# Define function 
f(x, y) = x^2 + y^2

# Plot the result 
fig, ax, plt = tridelaunayplotf(tridln, f, vmarkersize=10)
wireframe!(tomesh(tridln.delaunay), linewidth=3)
display(fig)

msh = tomesh(tridln.delaunay, f) 
pts = [Point2(val[1], val[2]) for val in msh.position]
msh2 = GeometryBasics.Mesh(pts, faces(msh))
mesh(msh) 
mesh!(msh2)
