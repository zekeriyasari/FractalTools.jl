# This file contains some of the utility functions.

using QuadGK

unit_step(t) = t == 0 ? 1 : (1/2*sign(t) + 1/2)

unit_pulse(t; tl=0, th=1) = (unit_step(t - tl) - unit_step(t - th))

fnorm1(f::Function, t0, tf) = quadgk(x -> abs(f(x)), t0, tf)[1]

fnormInf(f::Function, t0, tf) = maximum(f.(collect(range(t0, tf, length=1000))))

function fnormp(f::Function, t0, tf, p=1)
    if p == 1
        return fnorm1(f, t0, tf)
    elseif p == Inf
        return fnormInf(f, t0, tf)
    else
        return quadgk(x -> abs(f(x))^p, t0, tf)[1]^(1/p)
    end
end

ndgrid(v::AbstractVector) = copy(v)

function ndgrid(v1::AbstractVector, v2::AbstractVector)
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repeat(v1, 1, n), repeat(v2, m, 1))
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end

function ndgrid(vs::AbstractVector{T}...) where T
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Vector{T}(undef, sz), n)
    s = 1
    for i=1:n
        a = out[i]
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end

meshgrid(v::AbstractVector) = meshgrid(v, v)

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where T
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repeat(vx, m, 1), repeat(vy, 1, n))
end

function meshgrid(vx::AbstractVector, vy::AbstractVector, vz::AbstractVector)
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end
