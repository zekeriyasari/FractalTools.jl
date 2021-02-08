using Makie 

# Construct a figure window 
fig = Figure() 

# Construct two menus 
menu = Menu(fig, options = [:viridis, :heat, :blues])
funcs = [sqrt, x -> x^2, sin, cos] 
menu2 = Menu(fig, options=zip(["Square Root", "Square", "Sine", "Cosine"], funcs))

# Insert menus 
fig[1, 1] = vgrid!(
    Label(fig, "Colormap", width=nothing), 
    menu, 
    Label(fig, "Function", width=nothing), 
    menu2; 
    tellheight=false, width=200
) 

# Insert a an axis 
ax = Axis(fig[1, 2]) 
func = Node{Any}(funcs[1])
ys = @lift($func.(0 : 0.3 : 10))
scat = scatter!(ax, ys, color=ys) 

# Insert a color bar 
cb = Colorbar(fig[1, 3], scat, width=30)

# Bind interaction to widgets
on(menu.selection) do choice 
    scat.colormap = choice
end 

on(menu2.selection) do choice 
    func[] = choice
    autolimits!(ax)
end 

menu2.is_open = true

fig 

