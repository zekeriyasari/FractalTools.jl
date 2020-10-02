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
@bind d html"<input type=range min=0 max=1 step=1e-1>"

# ╔═╡ 2caf47f8-047f-11eb-1f15-db8e661dc9ce
@bind niter html"<input type=range min=1 max=100 step=1>"

# ╔═╡ 69dbf05e-047f-11eb-3d6e-892350bec24f
niter, d

# ╔═╡ 402f7b62-01b5-11eb-2789-7bdb04c0d694
begin
	# Consrucct data 
	xi, dx, xf = 0, 0.1, 1 
	x = collect(xi : dx : xf) 
	y = sin.(2π * x) 
	z = cos.(2π * x)
	dvar = d * ones(length(x) - 1)
	h = d * ones(length(x) - 1)
	l = d * ones(length(x) - 1)
	m = d * ones(length(x) - 1)
	
	# Interplate
	interp = fif(x, y, dvar, niter=niter) 
	
	# Plots 
	plt = plot()
	xd = collect(xi : 0.1dx : xf) 
	yd = interp.(xd)
	plot!(xd, yd)
	scatter!(x, y, marker=(:circle, 3))
end

# ╔═╡ 927c79ea-0482-11eb-2562-d1934e61f65e
@bind d_hidden html"<input type=range min=0 max=1 step=1e-1>"

# ╔═╡ 99eeb9cc-0482-11eb-3df1-b3a85faad5a6
@bind niter_hidden html"<input type=range min=1 max=100 step=1>"

# ╔═╡ 1078363c-0482-11eb-3b0a-e10a7b7dbd92
begin
	# Consrucct data 
	xi_hidden, dx_hidden, xf_hidden = 0, 0.1, 1 
	x_hidden = collect(xi_hidden : dx_hidden : xf_hidden) 
	y_hidden = sin.(2π * x_hidden) 
	z_hidden = cos.(2π * x_hidden)
	dvar_hidden = d * ones(length(x_hidden) - 1)
	h_hidden = d_hidden * ones(length(x_hidden) - 1)
	l_hidden = d_hidden * ones(length(x_hidden) - 1)
	m_hidden = d_hidden * ones(length(x_hidden) - 1)
	
	# Interplate
	interp_hidden = hiddenfif(x_hidden, y_hidden, z_hidden, dvar_hidden, h_hidden, l_hidden, m_hidden,niter=niter_hidden) 
	
	# Plots 
	plt_hidden = plot()
	xd_hidden = collect(xi_hidden : 0.1dx_hidden : xf_hidden) 
	yd_hidden = interp_hidden.(xd_hidden)
	plot!(xd_hidden, getindex.(yd_hidden, 1))
	scatter!(x_hidden, y_hidden, marker=(:circle, 3))
end

# ╔═╡ Cell order:
# ╠═2069c2c0-01b4-11eb-3cc9-a7c60467a92a
# ╠═24b8a340-01b5-11eb-0551-318932fe4a10
# ╠═2caf47f8-047f-11eb-1f15-db8e661dc9ce
# ╠═69dbf05e-047f-11eb-3d6e-892350bec24f
# ╠═402f7b62-01b5-11eb-2789-7bdb04c0d694
# ╠═927c79ea-0482-11eb-2562-d1934e61f65e
# ╠═99eeb9cc-0482-11eb-3df1-b3a85faad5a6
# ╠═1078363c-0482-11eb-3b0a-e10a7b7dbd92
