# This file includes the iterated function system (IFS) toos. 

export IFS, DetAlg, RandAlg, Attractor, attractor, detalg, randalg

""" 
    $(TYPEDEF) 

Iterated fucntion sytem (IFS) 

# Fields 

    $(TYPEDFIELDS)
"""
struct IFS{T1<:AbstractVector{<:Transformation}, T2<:AbstractVector{<:Real}}
    "Vector of transformations of IFS"
    ws::T1 
    "Vector of probabilities of IFS"
    probs::T2
    function IFS(ws::T1, probs::T2) where {T1, T2} 
        # Note: For the floating point numbers, aproximation(≈), instead of exact equal (==), should be considered
        sum(probs) ≈ 1 || throw(ArgumentError("Sum of probabilities must be 1."))
        new{T1, T2}(ws, probs) 
    end
end
IFS(ws) = (n = length(ws); IFS(ws, 1  / n * ones(n))) 

"""
    dimension(ifs::IFS)

Returns dimension of `ifs`.
"""
dimension(ifs::IFS) = dimension(ifs.ws[1])

"""
    contfactor(ifs::IFS)

Returns the contraction factor of `IFS`.
"""
contfactor(ifs::IFS) = maximum(contfactor.(ifs.ws))

"""
    $TYPEDEF

A type signifying that deterministic algorithm is used when calculating the attractor of and IFS.
"""
struct DetAlg end

"""
    $TYPEDEF

A type signifying that random algorithm is used when calculating the attractor of and IFS.
"""
struct RandAlg end

"""
    $TYPEDEF

Attractor of `IFS` type 

# Fields

    $TYPEDFIELDS
"""
struct Attractor{T, S, R}
    "IFS of Attractor"
    ifs::T
    "Type of algorithm to be used to compute attractor(Options are DetAlg and RandAlg"
    alg::S
    "Initial set of attractor"
    initset::R
    "Set of the attractor"
    set::R
    "Number of iterations"
    numiter::Int
    "Sequential or parallel"
    parallel::Bool
end

"""
    $SIGNATURES

Computes the attractor of `ifs`. If `alg` is of type `DetAlg`, the deterministic algorithm is used. If `alg` is of type
`RandAlg`, random algorithm is used. `kwargs` may include

* `numiter::Int` : Number of iterations to used to calcuate the attractor (defaults to 10)

* `numtransient::Int` : Number of transient iterations to used to calcuate a transient set. When the transient set is
  constructed, the computation of attractor is continued with distributed computation if `alg` is `RandAlg` and `parallel` is
  `true`. (defaults to 10)

* `parallel::Bool`: If  `true`, the attractor is computed using distrbuted computation. (defaults to false)

* `placedependent::Bool` : If `true`, place dependent attractor is computed if α and β is given accordingly. (default to
  false)

* `α::AbstractVector` : Place-dependent probility coefficient(defaults to nothing)

* `β::AbstractVector` : Place-dependent probility coefficient. (default to nothing)
"""
attractor(ifs, initset; alg=DetAlg(), kwargs...) = typeof(alg) == DetAlg ? 
                                                   detalg(ifs,initset; kwargs...) : 
                                                   randalg(ifs,initset; kwargs...)

"""
    $SIGNATURES

Computes the attractor of `ifs` with deterministic algorithm.`numiter` is number of iterations. (Defaults to 10). If
`parallel` is true, attractor is computed via parallel computation.
"""
function detalg(ifs, initset; numiter=10, parallel=false)
    copiedset = copy(initset)
    set = parallel ? 
          detalg_parallel(ifs.ws, copiedset, numiter) : 
          detalg_sequential(ifs.ws, copiedset, numiter)
    Attractor(ifs, DetAlg(), initset, set, numiter, parallel)
end

# Computes the attractor of an ifs via deterministic algorithm sequentially. 
function detalg_sequential(ws, set, numiter)
    for i in 1 : numiter
        set = vcat(map(w -> w.(set), ws)...)
    end
    set
end

# Computes the attractor of an ifs via deterministic algorithm in parallel. 
function detalg_parallel(ws, set, numiter)
    loadprocs()
    for i in 1 : numiter
        set = vcat(map(w -> pmap(w, set), ws)...)
    end
    set
end

"""
    $SIGNATURES

Computes the attractor of `ifs` with random algorithm.`numiter` is number of iterations. (Defaults to 100). `numtransient` is
the number of transient iterations. If `parallel` is true, attractor is computed via parallel computation. If
`placedependent` is true, the probabilties of the ifs are dependent on the coordinates `x`. This dependency `p(x)` is given
via the parameters `α` and `β` where p(x) = α x + β.
"""
function randalg(ifs, initset; numiter=100, numtransient=10, parallel=false, placedependent=false, α=nothing, β=nothing)
    ws = ifs.ws
    probs = ifs.probs
    if parallel
        if placedependent
            transient = randalg_sequential_pd(ws, copy(initset), numtransient, probs, α, β)
            set = randalg_parallel_pd(ws, transient, numiter, probs, α, β)
        else
            transient = randalg_sequential(ws, copy(initset), numtransient, probs)
            set = randalg_parallel(ws, transient, numiter, probs)
        end
    else
        if placedependent
            set = randalg_sequential_pd(ws, copy(initset), numiter, probs, α,β)
        else
            set = randalg_sequential(ws, copy(initset), numiter, probs)
        end
    end
    Attractor(ifs, RandAlg(), initset, set, numiter, parallel)
end

# Computes the attractor of an ifs via random algorithm sequentially. 
function randalg_sequential(ws, set, numiter, probs)
    weights = Weights(probs)
    xi = set[end]
    for i = 1 : numiter
        trfmi = sample(ws, weights)
        xi = trfmi(xi)
        push!(set, xi)
    end
    set
end


# Computes the attractor of an ifs via random algorithm sequentially with placedependent probabilties.
function randalg_sequential_pd(ws, set, numiter, probs, α, β)
    xi = set[end]
    for i = 1 : numiter
        trfmi = sample(ws, Weights(probs))
        xi = trfmi(xi)
        probs = α * xi + β
        push!(set, xi)
    end
    set
end

# Computes the attractor of an ifs via random algorithm in parallel. 
function randalg_parallel(ws, set, numiter, probs)
    weights = Weights(probs)
    loadprocs()
    vcat(pmap(process_chunk, [(ws, set, floor(Int, numiter / nworkers()), weights) for i =  1 : nworkers()])...)
end

# Computes the attractor of an ifs via random algorithm in parallel with placedependent probabilties.
function randalg_parallel_pd(ws, set, numiter, probs, α, β)
    loadprocs()
    vcat(pmap(process_chunk_pd, [(ws, set, floor(Int, numiter / nworkers()), probs, α, β) for i =  1 : nworkers()])...)
end

# `process_chunk` is the worker function that is used in all processes(both in master and worker process). when calculating
# the attractor if alg is `RandAlg` and `parallel` is true.
function process_chunk(ws_set_niter_weights)
    ws, set, niter, weights = ws_set_niter_weights
    xi = set[end]
    for i = 1 : niter
        wi = sample(ws, weights)
        xi = wi(xi)
        push!(set, xi)
    end
    set
end

# `process_chunk` is the worker function that is used in all processes(both in master and worker process). when calculating
# the attractor if alg is `RandAlg` and `parallel` is true with placedependent probabilties.
function process_chunk_pd(ws_set_niter_probs_alpha_beta)
    ws, set, niter, probs, α, β = ws_set_niter_probs_alpha_beta
    xi = set[end]
    for i = 1 : niter
        wi = sample(ws, Weights(probs))
        xi = wi(xi)
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


