export project, triangulate, findouttriangle, interpolate, gettransform 

"""
    $SIGNATURES

Returns a fractal surface interpolation function that interpolates `pnts3d`.
"""
function interpolate(pnts3d::AbstractVector{<:AbstractPoint}; α=0.01, f0 = (x, y) -> 0., maxiters=15, gettransforms::Bool=false)
    # Fint convex hull, i.e, the boundary triangle 
    outtrig = findouttriangle(pnts3d)

    # Consruct a 3d mesh 
    tess, msh3 = triangulate(pnts3d)
	domains = project(msh3)

    # Compute transforms 
    transforms = map(intrig -> gettransform(outtrig, intrig, α), msh3)
	ifs = IFS(transforms)

    # Define mapping 
    # footprints = Point2{Float64}[]      # Uncomment to make a container for the footprints of any point 
    function mapping(f)
		function fnext(xd, yd)
            pnt = Point(xd, yd)
            # push!(footprints, pnt)
			idx = tess.find_simplex(pnt)[1] + 1  
            idx == 0 && return NaN 
            transform = transforms[idx]
            A, b = transform.A, transform.b
			linv = A[1:2, 1:2] \ (pnt - b[1:2])
			A[3, 1 : 2] ⋅ linv + α * f(linv[1], linv[2]) + b[3]
        end
	end 

	# gettransforms ? (interpolant, transforms, footprints) : (interpolant, footprints)
    # Main iteration of the fractal interpolation.
    itp = ∘((mapping for i in 1 : maxiters)...)(f0)
	gettransforms ? (itp, transforms) : itp
	Interpolant(ifs, domains, itp)
end

# interpolate(pnts3d::AbstractVector; kwargs...) = interpolate(map(item -> Point(item...), pnts3d); kwargs...)

"""
    $SIGNATURES

Returns the two-dimensional projections of the three-dimensional points `pnts3d`
"""
project(pnts3d::AbstractVector{<:Point}) = [Point(pnt[1], pnt[2]) for pnt in pnts3d]

"""
    $SIGNATURES

Returns a two dimensional mesh whose points are the projections of `msh3`. 
"""
project(msh3::GeometryBasics.Mesh) = GeometryBasics.Mesh(project(msh3.position), copy(faces(msh3)))

"""
    $SIGNATURES

Returns a tuple of a Delaunay tessellation and a correspoding three dimensional mesh. 
"""
function triangulate(pnts3d) 
    pnts2d = project(pnts3d)
    tess = spt.Delaunay(pnts2d)
    trifaces = [TriangleFace(val...) for val in eachrow(tess.simplices .+ 1)]
    tess, GeometryBasics.Mesh(pnts3d, trifaces)
end 

"""
    $SIGNATURES 

Returns the bounding triangle of the points `pnts3d`.
"""
function findouttriangle(pnts3d)
    pnts2d = project(pnts3d)
    # TODO: Call c code directly
    # TODO: Check convecity and divide into convex dokmains.
    hull = spt.ConvexHull(collect(hcat(collect.(pnts2d)...)') )
    if length(hull.vertices) == 3
        # If the convex hull is a triangle, just return the triangle 
        Triangle(Point.(pnts3d[hull.vertices .+ 1])...)
    else 
        # If the the convex hull is not a triangle but a polygon, construct the 
        # boundary polygon, construct a mesh from the polygon and return the 
        # maximum triangle with the maximum area.
        polygon = Ngon(
            SVector{length(hull.vertices)}(Point.(pnts3d[hull.vertices .+ 1]))
        )
        msh = GeometryBasics.mesh(polygon)
        idx = argmax([area(trig.points) for trig in msh])
        msh[idx]
    end 
end

"""
    $SIGNATURES

Returns a named tuple of `A` and` `b` such that `L(x) = A * x + b` maps `outtrig` to `intrig`. 
"""
function gettransform(outtrig, intrig, α::Real=1.) 
    outmat = collect(hcat(coordinates(outtrig)...)')
    inmat = collect(hcat(coordinates(intrig)...)')
    inmat[:, end] -= α * outmat[:, end]
    outmat[:, end] = ones(3)
    sol = outmat \ inmat
	A = [  sol[1, 1]   sol[2, 1]   0;
            sol[1, 2]   sol[2, 2]   0; 
            sol[1, 3]   sol[2, 3]   α           ]
    b = [   sol[3, 1],  sol[3, 2],  sol[3, 3]   ]
	Transformation(A, b)
end 




# # This file includes two dimensional fractal interpolation functions.

# export fis

# """
#     fis(mesh, z, α=0.025, func0=(xi,yi)->0., num_iter=15, get_coeffs=false)

# Fractal interpolaton surface for the interpolaton domain given with `mesh` and
# associated `z` with the vertical scling factors `α`, initial function `func0`,
# number of iterations `num_iter`. If `get_coeffs` is `true`, coefficients of the
# IFS corresponding to dataset `mesh` and `z` are returned.
# """
# function fis(mesh::PyCall.PyObject, z::AbstractVector{<:Real}; α::Real=0.025, func0::Function=(xi,yi)->0.,
# 	num_iter::Int=15, get_coeffs::Bool=true, region_type::String="triangular")
# 	# Compute coefficients for each of the triangle in the mesh.
# 	coeffs = coefficients(mesh, z, α, region_type=region_type)

# 	# Get finder of the mesh
# 	get_finder = mesh.get_trifinder
# 	finder = get_finder()

#     # Construct nearest-neigbor tree of the mesh.
# 	tree = KDTree(collect(transpose([mesh.x mesh.y])))

# 	# Define template function for the main iteration.
# 	function mapping(func::Function)
# 		function fnext(xd::Real, yd::Real)
# 			tri_idx = finder(xd, yd)[1] + 1
#             # If (xd, yd) cannot be found, find the nearest triangle
#             # and project (xd, yd) back onto the triangle
# 			if tri_idx == 0
# 				nn_idx, nn_dist = knn(tree, [xd; yd], 1)
#                 tri_idx = finder(mesh.x[nn_idx], mesh.y[nn_idx])[1] .+ 1
# 				nn_tri_idx = mesh.triangles[tri_idx, :] .+ 1
# 				triang = collect(transpose([mesh.x[nn_tri_idx] mesh.y[nn_tri_idx]]))
# 				xd, yd = proj_onto_triangle([xd; yd], triang)
# 			end
#             A = coeffs[:, :, tri_idx]
# 			linv = collect(transpose(A[1:2, 1:2]) \ ([xd; yd] - A[3, 1:2]))
# 			return dot(linv, A[1:2, 3]) + α * func(linv[1], linv[2]) + A[3, 3]
# 		end
# 		return fnext
# 	end # function

# 	# Main iteration of the fractal mesh interpolation.
# 	interpolant = func0
# 	for i = 1 : num_iter
# 		interpolant = mapping(interpolant)
# 	end
# 	get_coeffs && return interpolant, coeffs
#     return interpolant
# end


# """
#     fis(mesh, z, α=0.1, func0=(xi,yi)->0., num_iter=15, get_coeffs=false)

# Fractal interpolaton surface for the interpolaton domain given with `mesh` and
# associated `z` with the vertical scling factors `α`, initial function `func0`,
# number of iterations `num_iter`. If `get_coeffs` is `true`, coefficients of the
# IFS orresponding to dataset `mesh` and `z` are returned.
# """
# function fis(regions::AbstractVector{PyCall.PyObject}, z::AbstractVector{<:AbstractVector}; 
# 	α::Real=0.1, func0::Function=(xi,yi)->0., num_iter::Int=10, get_coeffs::Bool=false)
# 	# Find piecewise interpolation functions for each region.
# 	num_regions = length(regions)
# 	funcs = Vector{Function}(undef, num_regions)  # Piecewise interpolation function.
# 	if get_coeffs
# 		coeffs = Vector{Array{Float64, 3}}(undef, num_regions)
# 		for k = 1 : num_regions
#             fs, cs = fis(regions[k], z[k], α=α, func0=func0, 
# 				num_iter=num_iter, get_coeffs=true)
# 			funcs[k] = fs
# 			coeffs[k] = cs
# 		end
# 		return (xd, yd)->funcs[find_region(regions, xd, yd)][1](xd, yd), coeffs
# 	end

# 	for k = 1 : num_regions
# 		funcs[k] = fis(regions[k], z[k], α=α, func0=func0, num_iter=num_iter, get_coeffs=false)
# 	end
# 	return (xd, yd) -> funcs[find_region(regions, xd, yd)][1](xd, yd)
# end

# """
#     coefficients(mesh, z, α, region_type)

# Computes coefficients of the iterated function system given the two dimensional
# interpolation domain `mesh` and corresponding `z` with the vertical scaling
# factor `α`.
# """
# function coefficients(mesh::PyCall.PyObject, z::AbstractVector{<:Real}, α::Real; 
# 	region_type::String="triangular", get_apex_coordinates::Bool=false)
#     # Note: For triagular regions, the apex points of the triangular region must be found to
#     # construct P matrix and zi vector
#     # For polygonal regions, the apex points of the maxiumum area triangle must be found to
#     # construct P matrix and zi vector.

#     if startswith(region_type, "triangular")
# 		j, apex_coords = find_apex(mesh, coordinates=true)
#         P = [mesh.x[j] mesh.y[j] ones(3, 1)]
#         zj = z[j]
#     elseif startswith(region_type, "polygonal")
#         # Find P matrix and its inverse.
#         _, apex_coords = find_apex(mesh, coordinates=true)
#         area, a, b, c = maximum_area_triangle(apex_coords)
#         acoord = apex_coords[a.idx, :]
#         bcoord = apex_coords[b.idx, :]
#         ccoord = apex_coords[c.idx, :]
#         P = [hcat(acoord, bcoord, ccoord)' ones(3, 1)]
#         zi = [z[a.idx], z[b.idx], z[c.idx]]
#     else
#         error("Unknown region type. Available region types are `triangular` and `polygonal`")
#     end

# 	# Compute coefficients for each triangle.
# 	num_triangles = size(mesh.triangles, 1)
# 	coeffs = zeros(3, 3, num_triangles)
# 	for i = 1 : num_triangles
# 		k = mesh.triangles[i, :] .+ 1  # Node indices
# 	    xi, yi, zi = mesh.x[k], mesh.y[k], z[k]
# 		coeffs[:, :, i] = P \ [xi yi (zi - α * zj)]
# 	end 
# 	if get_apex_coordinates
# 		return coeffs, apex_coords
# 	end
# 	return coeffs
# end 

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
