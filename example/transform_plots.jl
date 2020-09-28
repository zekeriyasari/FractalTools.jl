using FractalTools
using Plots 

# Define function 
f(x) = 4 - (x - 2)^2 + 1

# Interpolation points 
xi, dx, xf = 0, 1, 4 
x = collect(xi : dx : xf) 
y = f.(x) 
d = 0.5 * ones(length(x) - 1)

# Construct interpolation 
interp = fif(x, y, d)
ws = interp.ifs.ws

# End points 
x = [0, 4, 4, 0, 0]
y = [0, 0, 4, 4, 0]
p = [[x[i], y[i]] for i in 1 : length(x)]

function transform_plot()
    plt = plot(getindex.(p, 1), getindex.(p,2), marker=(:circle, 3))
    xt = collect(xi : 0.1dx : xf)
    plot!(xt, f.(xt))
    for wt in ws 
        wp = []
        for (xi, yi) in zip(x, y)
            push!(wp, wt([xi, yi]))
        end 
        plot!(getindex.(wp, 1), getindex.(wp,2), marker=(:circle, 3), subplot=1)
    end
    plt
end 

transform_plot()
