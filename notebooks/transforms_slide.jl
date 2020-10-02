### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ e6d6879e-048f-11eb-1328-2d8ae06f1f1f
# Include package. 
begin 
	using Pkg
	Pkg.activate(joinpath(Pkg.envdir(), "dev-env"))
	using Plots 
	using FractalTools 
end 

# ╔═╡ 2df6b090-0490-11eb-1093-ff3605ecef25
# Define function 
f(x) = 4 - (x - 2)^2 + 1

# ╔═╡ 3f42bc90-0490-11eb-3570-3bafd44e6d9d
# Interpolation points 
begin 
	xi, dx, xf = 0, 1, 4 
	x = collect(xi : dx : xf) 
	y = f.(x) 
	d = 0.5 * ones(length(x) - 1)	
end

# ╔═╡ 4abce91a-0490-11eb-0841-9d7d87776e2c
# Construct interpolation 
begin 
	interp = fif(x, y, d)
	ws = interp.ifs.ws
end

# ╔═╡ 595c870a-0490-11eb-12fa-09f822a8d78d

function transform_plot(niters)
	xc = [0, 4, 4, 0, 0]
	yc = [0, 0, 4, 4, 0]
	pc = [xc yc]
	xt = collect(xi : 0.1dx : xf)
	plt = plot(xt, f.(xt))
	for niter in 1 : niters
		pc = vcat([mapslices(wt, pc, dims=2) for wt in ws]...)
		plot!(pc[:, 1], pc[:, 2])
	end
	plt
end


# ╔═╡ cb8afd34-0490-11eb-2d48-b3718511cf5d
@bind niters html"<input type=range min=1 max=10 step=1>"

# ╔═╡ f7275e38-0490-11eb-333d-5be95205a2f4
niters

# ╔═╡ ce9dd820-0490-11eb-22da-f9ba42b2c68a
transform_plot(niters)

# ╔═╡ Cell order:
# ╠═e6d6879e-048f-11eb-1328-2d8ae06f1f1f
# ╠═2df6b090-0490-11eb-1093-ff3605ecef25
# ╠═3f42bc90-0490-11eb-3570-3bafd44e6d9d
# ╠═4abce91a-0490-11eb-0841-9d7d87776e2c
# ╠═595c870a-0490-11eb-12fa-09f822a8d78d
# ╠═cb8afd34-0490-11eb-2d48-b3718511cf5d
# ╠═f7275e38-0490-11eb-333d-5be95205a2f4
# ╠═ce9dd820-0490-11eb-22da-f9ba42b2c68a
