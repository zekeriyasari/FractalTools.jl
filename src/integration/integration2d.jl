# This file includes two dimensional integration.


"""
    integrate(mesh::PyCall.PyObject, z::AbstractVector; α::Real=0.001, region_type::String="triangular")

I = ∫f(x,y)dxdy

Numerically integrates interpolation function for interpolation points (`mesh`, `zi`)
where `mesh` is domain, and `zi` is the corresponding point in the
range. `α` is vertical scaling factor between (0., 1.). `region_type` determines the
geometry of the regions. `triangular` and `polygonal` are the available geometrical regions.
"""
function integrate(mesh::PyCall.PyObject, z::AbstractVector; α::Real=0.001, region_type::String="triangular")
    # Compute transformation coefficients
    coeffs, apex_coords = coefficients(mesh, z, α, region_type=region_type, get_apex_coordinates=true)
    α5 = coeffs[1, 3, :]
    α6 = coeffs[2, 3, :]
    β3 = coeffs[3, 3, :]

    # Jacobien of determinant of change of variables transformation. The change of variable
    # transformation is to transform smaller integration regions to maximum area triangular
    # region.
    ΔL = abs.(det.([coeffs[1:2, 1:2, i] for i = 1 : size(coeffs, 3)]))

    # `domain_vec` and `domain_mat` is the vector and matrix, respectively, of the
    # transformation U to transform each smaller triangle to unit triangle.
    # U(x, y) = [a11 a12; a12 a22 ][x; y] + [b1; b2] = A x + b
    xt1, xt2, xt3 = apex_coords[:, 1]
    yt1, yt2, yt3 = apex_coords[:, 2]
    @info (xt1, xt2, xt3)
    @info (yt1, yt2, yt3)
    # domain_vec = [0, 0, 1, 0, 0, 1]
    # domain_mat = [x1 y1 0 0 1 0;
    #               0 0 x1 y1 0 1;
    #               x2 y2 0 0 1 0;
    #               0 0 x2 y2 0 1;
    #               x3 y3 0 0 1 0;
    #               0 0 x3 y3 0 1]
    # @info domain_mat
    domain_vec = [xt1, yt1, yt2, xt3, xt2, yt3]
    domain_mat = [0 0 1 0 0 0;
                  0 0 0 0 0 1;
                  0 0 0 1 0 1;
                  0 1 1 0 0 0;
                  1 0 1 0 0 0;
                  0 0 0 0 1 1]
    A = domain_mat \ domain_vec
    ΔU = abs.(det([A[1] A[2];
                   A[3] A[4]]))
    k11 = A[1] + A[2] + 3*A[5]
    k12 = A[3] + A[4] + 3*A[6]
    K1 = ΔU/6 * sum(ΔL .* (α5 * k11 + α6 * k12 + 3*β3))
    K2 = sum(ΔL * α)
    I = K1 / (1 - K2)  # Integration value
    return I
end
