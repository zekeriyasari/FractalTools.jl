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

# ╔═╡ 2069c2c0-01b4-11eb-3cc9-a7c60467a92a
# Load packages 
begin 
	using Pkg 
	Pkg.activate(mktempdir())
	Pkg.add("Plots")
	Pkg.develop(path="/home/sari/.julia/dev/FractalTools/")
	using Plots
	plotly()
	using FractalTools
end 

# ╔═╡ 24b8a340-01b5-11eb-0551-318932fe4a10
@bind d html"<input type=range min=0 max=1 step=1e-2>"

# ╔═╡ 402f7b62-01b5-11eb-2789-7bdb04c0d694
begin
	# Consrucct data 
	xi, dx, xf = 0, 0.1, 1 
	x = collect(xi : dx : xf) 
	y = sin.(2π * x) 
	dvar = d * ones(length(x) - 1)
	
	# Interplate
	interp = fif(x, y, dvar) 
	
	# Plots 
	plt = plot()
	xd = collect(xi : 0.1dx : xf) 
	yd = interp.(xd)
	plot!(xd, yd)
	scatter!(x, y, marker=(:circle, 3))
end

# ╔═╡ Cell order:
# ╠═2069c2c0-01b4-11eb-3cc9-a7c60467a92a
# ╠═24b8a340-01b5-11eb-0551-318932fe4a10
# ╠═402f7b62-01b5-11eb-2789-7bdb04c0d694
