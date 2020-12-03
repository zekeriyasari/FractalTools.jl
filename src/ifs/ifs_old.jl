# This file includes IFS type of FractalTools

"""
    IFS(trfms, probs, set)

Constructs an IFS. `trfms` is the affine transformations of IFS. `probs` is the probabilities of `trfms`. `set` is the
initial set in the metric space of `IFS`.

    IFS(trfms)

Constructs and `IFS` from `trfms` with equal probabilities and a random initial set.
"""
mutable struct IFS{T<:AbstractAffineTransform}
    trfms::Vector{T}
    probs::Vector{Float64}
    set::Vector{Vector}
    function IFS(trfms::AbstractVector{T}, probs=ones(length(trfms))/length(trfms),
        set=[rand(dimension(trfms[1]))]) where T <: AbstractAffineTransform
        sum(probs) > 1 && error("Sum of probs must not be greater than 1.")
        new{eltype(trfms)}(trfms, probs, set)
    end
end

show(io::IO, ifs::IFS) = print(io, "IFS(numtrfms:$(length(ifs.trfms)), probs:$(ifs.probs), set:$(ifs.set))")

"""
    dimension(ifs::IFS)

Returns dimension of `ifs`.
"""
dimension(ifs::IFS) = dimension(ifs.trfms[1])

"""
    contfactor(ifs::IFS)

Returns the contraction factor of `IFS`.
"""
contfactor(ifs::IFS) = maximum(contfactor.(ifs.trfms))

"""
    DetAlg

A type signifying that deterministic algorithm is used when calculating the attractor of and IFS.
See also: [`attractor`](@ref)
"""
struct DetAlg end

"""
    RandAlg

A type signifying that random algorithm is used when calculating the attractor of and IFS.
See also: [`attractor`](@ref)
"""
struct RandAlg end

struct Attractor{T, S, R}
    ifs::T
    alg::S
    set::R
end

"""
    attractor(ifs, alg=DetAlg(), kwargs...)

Computes the attractor of `ifs`. If `alg` is `DetAlg`, the deterministic algorithm is used. If `alg` is `RandAlg`,
random algorithm is used. kwargs may include

* `numiter::Int` : Number of iterations to used to calcuate the attractor (defaults to 10)
* `numtransient::Int` : Number of transient iterations to used to calcuate a transient set. When the transient set is
    constructed, the computation of attractor is continued with distributed computation if `alg` is `RandAlg` and
    `parallel` is `true`. (defaults to 10)
* `parallel::Bool`: If  `true`, the attractor is computed using distrbuted computation. (defaults to false)
* `placedependent::Bool` : If `true`, place dependent attractor is computed if α and β is given accordingly. (default to false)
* `α::AbstractVector` : Place-dependent probility coefficient. (defaults to nothing)
* `β::AbstractVector` : Place-dependent probility coefficient. (default to nothing)
"""
function attractor(ifs; alg=DetAlg(), kwargs...)
    set = typeof(alg) == DetAlg ? detalg(ifs; kwargs...) : randalg(ifs; kwargs...)
    return Attractor(ifs, alg, set)
end

function detalg(ifs; numiter=10, parallel=false)
    trfms = ifs.trfms
    weights = Weights(ifs.probs)
    set = ifs.set
    if parallel
        detalg_parallel(trfms, set, numiter)
    else
        detalg_sequential(trfms, set, numiter)
    end
end

function detalg_sequential(trfms, set, numiter)
    for i in 1 : numiter
        set = vcat(map(trfm -> trfm.(set), trfms)...)
    end
    return set
end

function detalg_parallel(trfms, set, numiter)
    loadprocs()
    for i in 1 : numiter
        set = vcat(map(trfm -> pmap(trfm, set), trfms)...)
    end
    return set
end

function randalg(ifs; numiter=100, numtransient=10, parallel=false, placedependent=false, α=nothing, β=nothing)
    trfms = ifs.trfms
    probs = ifs.probs
    set = ifs.set
    if parallel
        if placedependent
            transient = randalg_sequential_pd(trfms, set, numtransient, probs, α, β)
            randalg_parallel_pd(trfms, transient, numiter, probs, α, β)
        else
            transient = randalg_sequential(trfms, set, numtransient, probs)
            randalg_parallel(trfms, transient, numiter, probs)
        end
    else
        if placedependent
            randalg_sequential_pd(trfms, set, numiter, probs, α,β)
        else
            randalg_sequential(trfms, set, numiter, probs)
        end
    end
end


function randalg_sequential(trfms, set, numiter, probs)
    weights = Weights(probs)
    xi = set[end]
    for i = 1 : numiter
        trfmi = sample(trfms, weights)
        xi = trfmi(xi)
        push!(set, xi)
    end
    return set
end

function randalg_sequential_pd(trfms, set, numiter, probs, α, β)
    xi = set[end]
    for i = 1 : numiter
        trfmi = sample(trfms, Weights(probs))
        xi = trfmi(xi)
        probs = α * xi + β
        push!(set, xi)
    end
    return set
end

function randalg_parallel(trfms, set, numiter, probs)
    weights = Weights(probs)
    loadprocs()
    vcat(pmap(process_chunk, [(trfms, set, floor(Int, numiter / nworkers()), weights) for i =  1 : nworkers()])...)
end

function randalg_parallel_pd(trfms, set, numiter, probs, α, β)
    loadprocs()
    vcat(pmap(process_chunk_pd, [(trfms, set, floor(Int, numiter / nworkers()), probs, α, β) for i =  1 : nworkers()])...)
end

# `process_chunk` is the worker fucntion that is used in all processes(both in master and worker process).
# when calculating the attractor if alg is `RandAlg` and `parallel` is true.
function process_chunk(trfms_set_niter_weights)
    trfms, set, niter, weights = trfms_set_niter_weights
    xi = set[end]
    for i = 1 : niter
        trfmi = sample(trfms, weights)
        xi = trfmi(xi)
        push!(set, xi)
    end
    set
end

function process_chunk_pd(trfms_set_niter_probs_alpha_beta)
    trfms, set, niter, probs, α, β = trfms_set_niter_probs_alpha_beta
    xi = set[end]
    for i = 1 : niter
        trfmi = sample(trfms, Weights(probs))
        xi = trfmi(xi)
        probs = α * xi + β
        push!(set, xi)
    end
    set
end

# Load worker processes and load FractalTools to those worker processes.
function loadprocs(numprocs=Base.Sys.CPU_THREADS - 1 - nprocs())
    addprocs(numprocs)
    @everywhere @eval using FractalTools
end

