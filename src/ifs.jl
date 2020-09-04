# This file includes the iterated function system (IFS) toos. 

export IFS 

""" 
    $(TYPEDEF) 

Iterated fucntion sytem (IFS) 

# Fields 

    $(TYPEDFIELDS)
"""
struct IFS{T1, T2}
    ws::T1 
    probs::T2
    IFS(ws::AbstractVector{T1}, probs::AbstractVector{T2}) where {T1<:Transformation, T2<:Real} = 
        new{typeof(ws), typeof(probs)}(ws, probs) 
end
IFS(ws) = (n = length(ws); IFS(ws, 1  / n * ones(n))) 