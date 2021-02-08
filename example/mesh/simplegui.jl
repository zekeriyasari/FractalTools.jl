# This file is a simple gui 
using Makie 

# Construct figure window
fig = Figure() 

# Construct a slider, a button and a plot axis
sl = Slider(fig[1, 1], range = 1 : 0.01 : 10)
ax = Axis(fig[2, 1])
btn = Button(fig[3, 1], tellwidth=false)

# Define a node and bind actions for interactivity
f = Node(1.)
t = collect(0. : 0.01 : 1.)
y = @lift(sin.(2π * $f * t))

# Define interactivity for slider 
on(sl.value) do val 
    f[] = val 
end 
# Define interactivty for button 
on(btn.clicks) do n 
    f[] += 1
end

foo = lines!(ax, t, y, linewidth = 1)

fig 