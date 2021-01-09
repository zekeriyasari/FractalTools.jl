# This file includes two dimensional fractal interpolation functions.

export hiddenfis

"""
    hiddenfis(mesh, z, α=0.025, func0=(xi,yi)->0., num_iter=15, get_coeffs=false)

Fractal interpolaton surface for the interpolaton domain given with `mesh` and
associated `z` with the vertical scling factors `α`, initial function `func0`,
number of iterations `num_iter`. If `get_coeffs` is `true`, coefficients of the
IFS corresponding to dataset `mesh` and `z` are returned.
"""
function hiddenfis(mesh::PyCall.PyObject, z::AbstractVector{<:Real}, t::AbstractVector{<:Real}; 
    α7::Real=0.025, α8::Real=0.025, α11::Real=0.025, α12::Real=0.025, func0::Function=(xi,yi) -> [0., 0.],
	num_iter::Int=15, get_coeffs::Bool=true, region_type::String="triangular")
	# Compute coefficients for each of the triangle in the mesh.
	coeffs = coefficients(mesh, z, t, α7, α8, α11, α12, region_type=region_type)

	# Get finder of the mesh
	get_finder = mesh.get_trifinder
	finder = get_finder()

    # Construct nearest-neigbor tree of the mesh.
	tree = KDTree(collect(transpose([mesh.x mesh.y])))

	# Define template function for the main iteration.
	function mapping(func::Function)
		function fnext(xd::Real, yd::Real)
			tri_idx = finder(xd, yd)[1] + 1
            # If (xd, yd) cannot be found, find the nearest triangle
            # and project (xd, yd) back onto the triangle
			if tri_idx == 0
				@info xd, yd
				# nn_idx, nn_dist = knn(tree, [xd; yd], 1)
                # tri_idx = finder(mesh.x[nn_idx], mesh.y[nn_idx])[1] .+ 1
				# nn_tri_idx = mesh.triangles[tri_idx, :] .+ 1
				# triang = collect(transpose([mesh.x[nn_tri_idx] mesh.y[nn_tri_idx]]))
				# xd, yd = proj_onto_triangle([xd; yd], triang)
				return [NaN, NaN]
			end
			A = coeffs[:, :, tri_idx]
			α = [α7 α8; α11 α12]; 
			linv = collect(transpose(A[1:2, 1:2]) \ ([xd; yd] - A[3, 1:2]))
			return transpose(A[1:2, 3:4]) * linv + α * func(linv[1], linv[2]) + A[3, 3:4]
		end
		return fnext
	end # function

	# Main iteration of the fractal mesh interpolation.
	interpolant = func0
	for i = 1 : num_iter
		interpolant = mapping(interpolant)
	end
	get_coeffs && return interpolant, coeffs
    return interpolant
end


"""
    hiddenfis(mesh, z, α=0.1, func0=(xi,yi)->0., num_iter=15, get_coeffs=false)

Fractal interpolaton surface for the interpolaton domain given with `mesh` and
associated `z` with the vertical scling factors `α`, initial function `func0`,
number of iterations `num_iter`. If `get_coeffs` is `true`, coefficients of the
IFS orresponding to dataset `mesh` and `z` are returned.
"""
function hiddenfis(regions::AbstractVector{PyCall.PyObject}, z::AbstractVector{<:AbstractVector}, 
	t::AbstractVector{<:AbstractVector}; 
	α7::Real=0.1, α8::Real=0.1, α11::Real=0.1, α12::Real=0.1 , func0::Function=(xi,yi)->[0., 0.], num_iter::Int=10, 
	get_coeffs::Bool=false)
	# Find piecewise interpolation functions for each region.
	num_regions = length(regions)
	funcs = Vector{Function}(undef, num_regions)  # Piecewise interpolation function.
	if get_coeffs
		coeffs = Vector{Array{Float64, 3}}(undef, num_regions)
		for k = 1 : num_regions
            fs, cs = hiddenfis(regions[k], z[k], t[k], α7=α7, α8=α8, α11=α11, α12=α12, func0=func0, 
				num_iter=num_iter, get_coeffs=true)
			funcs[k] = fs
			coeffs[k] = cs
		end
		return (xd, yd)->funcs[find_region(regions, xd, yd)][1](xd, yd), coeffs
	end

	for k = 1 : num_regions
		funcs[k] = hiddenfis(regions[k], z[k], t[k], α7=α7, α8=α8, α11=α11, α12=α12, func0=func0, num_iter=num_iter, get_coeffs=false)
	end
	return (xd, yd) -> funcs[find_region(regions, xd, yd)][1](xd, yd)
end

"""
    coefficients(mesh, z, α, region_type)

Computes coefficients of the iterated function system given the two dimensional
interpolation domain `mesh` and corresponding `z` with the vertical scaling
factor `α`.
"""
function coefficients(mesh::PyCall.PyObject, z::AbstractVector{<:Real}, t::AbstractVector{<:Real}, α7::Real, α8::Real, α11::Real, α12::Real; 
	region_type::String="triangular", get_apex_coordinates::Bool=false)
    # Note: For triagular regions, the apex points of the triangular region must be found to
    # construct P matrix and zi vector
    # For polygonal regions, the apex points of the maxiumum area triangle must be found to
    # construct P matrix and zi vector.

    if startswith(region_type, "triangular")
		j, apex_coords = find_apex(mesh, coordinates=true)
        P = [mesh.x[j] mesh.y[j] ones(3, 1)]
        zj = z[j]
        tj = t[j]
    elseif startswith(region_type, "polygonal")
        # Find P matrix and its inverse.
        _, apex_coords = find_apex(mesh, coordinates=true)
        area, a, b, c = maximum_area_triangle(apex_coords)
        acoord = apex_coords[a.idx, :]
        bcoord = apex_coords[b.idx, :]
        ccoord = apex_coords[c.idx, :]
        P = [hcat(acoord, bcoord, ccoord)' ones(3, 1)]
        zi = [z[a.idx], z[b.idx], z[c.idx]]
    else
        error("Unknown region type. Available region types are `triangular` and `polygonal`")
    end

	# Compute coefficients for each triangle.
	num_triangles = size(mesh.triangles, 1)
	coeffs = zeros(3, 4, num_triangles)
	for i = 1 : num_triangles
		k = mesh.triangles[i, :] .+ 1  # Node indices
	    xi, yi, zi, ti = mesh.x[k], mesh.y[k], z[k], t[k]
		coeffs[:, :, i] = P \ [xi yi (zi - α7 * zj - α8 * tj) (ti - α11 * zj - α12 * tj)]
	end 
	if get_apex_coordinates
		return coeffs, apex_coords
	end
	return coeffs
end 

# """
#     proj_onto_line(px, p1, p2, minimizer::Bool=true)

# Projects the point `px` onto the line passing the points `p1` ans `p2`.
# If `back_project` is `true`, minimizer of the disatance of the point to the line
# is returned instead of the projection point.
# """
# function proj_onto_line(px, p1,p2, minimizer=true)
# 	α = dot((p2 - p1), (px - p1)) / norm(p2 - p1)^2
# 	if minimizer
# 		α = min(1, max(α, 0))
# 	end
# 	return p1 + α * (p2 - p1)
# end

# """
#     proj_onto_triangle(px, triangle, all)

# Projects to the poits `px` on the nearest boundary of the triangle `triangle`.
# If `all` is `true`, then the projections of the point `px` is all three edges
# of the `triangle` is returned.
# """
# function proj_onto_triangle(px, triangle, all=false)
#     pα = hcat([proj_onto_line(px, triangle[:, idx[1]], triangle[:, idx[2]])  for idx in combinations(1:3, 2)]...)
# 	err = abs.(pα .- px)
# 	minval, minidx = findmin([norm(err[:, i]) for i = 1 : 3])
# 	if all
# 		return pα
# 	end
# 	return pα[:, minidx]
# end

# """
#     find_region(regions, px, py)

# Finds triangular region in `regions` to which the point `px` belongs.
# """
# function find_region(regions, px, py)
# 	num_regions = length(regions)
# 	mask = falses(num_regions)
# 	for (i, mesh) in enumerate(regions)
# 		get_finder = mesh.get_trifinder
# 		finder = get_finder()
# 		mask[i] = finder(px, py)[1] >= 0
# 	end

# 	# If (px, py) is on the boundary of the regions so that it could not be
# 	# be found, then find the region of its nearest neighbor.
# 	if all(mask .== false)
# 		x = vcat([region.x for region in regions]...)
# 		y = vcat([region.y for region in regions]...)
# 		tree = KDTree(collect(transpose([x y])))
# 		idx, dist = knn(tree, [px; py], 1)
# 		return find_region(regions, x[idx][1], y[idx][1])
# 	end
# 	return mask
# end 

# """
#     find_apex(mesh, coordinates=false)

# Finds the indices of the apex points of the `mesh`. If `coordinates` is true,
# the coordinates of the apex points are returned.
# """
# function find_apex(mesh; coordinates=false)
# 	convex_hull = spt.ConvexHull
# 	idx = convex_hull([mesh.x mesh.y]).vertices .+ 1  # Apex indices.
# 	if coordinates
# 		coords = [mesh.x[idx] mesh.y[idx]]  # Apex coordinates.
# 		return idx, coords
# 	end
# 	return idx
# end 

# """
#     contraction_factors(coeffs, α)

# Computes contraction factors of the affine transformations of the ifs given with
# `coeffs` and the verticaling factors `α`. Contraction factors are calculated
# to be the eigenvalues of the `coeffs`.
# """
# function contraction_factors(coeffs, α)
# 	num_factors = size(coeffs, 3)
#     factors = zeros(num_factors)
# 	for k = 1 : num_factors
# 		A = [coeffs[1, 1, k] coeffs[2, 1, k] 0.;
# 		 	 coeffs[1, 2, k] coeffs[2, 2, k] 0.;
# 			 coeffs[1, 3, k] coeffs[2, 3, k] α]
# 		# contraction_factors[k] = norm(A, p_norm)
# 		# contraction_factors[k] = abs(det(A))
# 		factors[k] = maximum(abs(eigvals(A)))  # Contraction factors.
# 	end
# 	return maximum(factors)
# end

# """
#     coefficients2equationsystem(coeffs, α)

# Returns the matrix, vector corresponding to the IFS given by `coeffs` and
# the vertical scaling factors `α`.
# """
# function coefficients2equationsystem(coeffs, α)
# 	A = [coeffs[1, 1] coeffs[2, 1] 0.;
# 		 coeffs[1, 2] coeffs[2, 2] 0.;
# 		 coeffs[1, 3] coeffs[2, 3] α]
# 	b = coeffs[3, :]
# 	return A, b
# end

# """
#     coefficients2transforms(coeffs, α)

# Converts to equation system given with `coeffs` and vertical scaling factors
# `α7` into transformations of the IFS.
# """
# function coefficients2transforms(coeffs::AbstractArray{<:Real, 2}, α)
# 	A, b = coefficients2equationsystem(coeffs, α)
# 	return x -> A * x + b
# end

# """
#     coefficients2transforms(coeffs, α)

# Converts to equation system given with `coeffs` and vertical scaling factors
# `α7` into transformations of the IFS.
# """
# function coefficients2transforms(coeffs::AbstractArray{<:Real, 3}, α)
# 	num_transforms = size(coeffs, 3)
# 	transforms = Array{Function}(num_transforms)
# 	for k = 1 : num_transforms
# 		transforms[k] = coefficients2transforms(coeffs[:, :, k], α)
# 	end
# 	return transforms
# end

# function find_vertex_neighbours(mesh, v)
#     x, y = v
#     set_mask = mesh.set_mask
#     get_finder = mesh.get_trifinder
#     finder = get_finder()
#     mask = falses(size(mesh.triangles, 1))  # Initially no mosk.
#     set_mask(mask)
#     # Get finder of the mesh
#     idxs = zeros(Int, 6)
#     for i = 1 : length(idxs)
#         tri_idx = finder(x, y)[1] + 1
#         if tri_idx != 0
#             idxs[i] = tri_idx
#             mask[tri_idx] .= true
#             set_mask(mask)
#         else
#             @warn "$tri_idx not found."
#             continue
#         end
#     end
#     idxs[idxs .!= 0]
# end

# function reset_mask(mesh, mask = nothing)
#     set_mask = mesh[:set_mask]
#     if mask === nothing
#         mask = falses(size(mesh[:triangles], 1))  # Initially no mosk
#     end
#     set_mask(mask)
# end
