# This file includes an example to reset the camera view. 

xs = LinRange(0, 10, 100)
ys = LinRange(0, 15, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]
fig, ax, plt = surface(xs, ys, zs)

# Lines below resets the camera view.
update_cam!(ax.scene, Vec3f0(3), Vec3f0(0))
center!(ax.scene)

