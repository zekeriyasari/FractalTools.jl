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

include("transformations.jl")
include("ifs.jl")
include("fif.jl")
include("hiddenfif.jl")
include("fis.jl")

function __init__()
    global spt = pyimport("scipy.spatial")
end

end # module
