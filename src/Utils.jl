## This file contains utility functions to work with polynomials
##
##

hstack(A) = reduce(hcat, A)
vstack(A) = reduce(vcat, A)

## - Polynomials are represented as vectors of coefficients.
## - Arrays of polynomials are represented as a single vector of the concatenated coefficients.
complexfy(v) = [x + y*im for (x,y) in zip(v[begin:2:end], v[(begin+1):2:end])]
splitDIM(curve) = collect.(Iterators.partition(curve, DIM))
curve_to_polys(curve) = splitDIM(complexfy(curve))
# curve_to_polys(curve::AbstractVector{<:Complex}) = curve_to_polys(curve)
evalcurve(curve::AbstractVector{<:Real}, coordinate::T) where {T<:Number} =
    evalpoly.([coordinate], curve_to_polys(curve))

randcurve(::Type{T}=Float64) where {T} = randn(T, N * (2 * DIM))

function set_common!(p_out, p)
    for i in eachindex(p)
        p_out[i] = p[i]
    end
    return p_out
end

_derivative(poly::AbstractVector) = poly[(begin+1):end].*(1:(length(poly)-1))
function _derivative(poly::AbstractVector, n)
    if n < 2
        return _derivative(poly)
    end
    return _derivative(_derivative(poly), n-1)
    # return _derivative(__derivative(v, n-1))
end
# _distancecurve(curve1, curve2) = maximum(LA.norm.(curve1 .- curve2); init = 0.)

## Some functions for projetive vectors
# vectors = [v1,...,vn] with length(vi) = length(vj) => length(vectors)
# So, v1,...,vn are the rows of a matrix, and we have less vectors than dimensions.
maxminors(vectors...) = zip(CC.combinations.(vectors, length(vectors))...)
isnull((px, py), (qx, qy); kwargs...) = isapprox(px * qy, py * qx; kwargs...)
equal_projective(P, Q; kwargs...) = all(minor -> isnull(minor...; kwargs...), maxminors(P, Q))
_dehomo(ppoint) = ppoint[begin:(end -1)]./last(ppoint)


# Fuzzy "uniquefy" an array.
intol(a, A; kwargs...) = any(x -> isapprox(x, a; kwargs...), A)
function intersectol(f::Function, A::AbstractArray{T}, B::AbstractArray{S}; kwarg...) where {T, S}
    out = T[]
    for a in A
        b = f(a)
        if intol(b, B; kwarg...)
            push!(out, a)
        end
    end
    return out
end
intersectol(A, B; kwarg...) = intersectol(identity, promote(A, B)...; kwarg...)

function uniquetol(f::Function, A::AbstractArray{T}, ::Type{S}=typeof(f(first(A))); kwargs...) where {T, S}
    out = T[]
    images = S[]
    for a in A
        b = f(a)
        if !(intol(b, images; kwargs...))
            push!(images, b)
            push!(out, a)
        end
    end
    return out
end
uniquetol(A; kwargs...) = uniquetol(identity, A, eltype(A); kwargs...)

real1(X...) = <((real∘first).(X)...) # real(first(x)) < real(first(y))
sort_byreal!(V) = sort!(V, lt = real1)
_flip(pair) = [last(pair), first(pair)]
