module FractalTools

import Base: show, display
import StatsBase: sample, Weights
using DocStringExtensions
using LinearAlgebra 
using Distributed
using Combinatorics
using PyCall
using HCubature
using LinearAlgebra
using NearestNeighbors
using GeometryBasics
using Clustering

include("mesh.jl")

include("ifs/transformations.jl")
include("ifs/ifs.jl")
include("ifs/prototypes.jl")
include("interpolation/fif.jl")
include("interpolation/hiddenfif.jl")
include("interpolation/fis.jl")
include("interpolation/hiddenfis.jl")
include("integration/integration1d.jl")
include("integration/integration2d.jl")

function __init__()
    global spt = pyimport_conda("scipy.spatial", "scipy")
end

end # module
