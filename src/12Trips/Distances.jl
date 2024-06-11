#####################
### Fubini Study distance in complex projective plane PC^2
######################
## Assume 'p' and 'q' have unitary length
## To bypass 'acos's safety error 'abs(x)>1'
_acos(x) = abs(x) <= 1 ? acos(x) : zero(x)
_FubiniStudy(z::T) where {T<:Number} = (_acos∘LA.norm)(z)
_FubiniStudy(p, q) = _FubiniStudy(LA.dot(p, q))

@doc raw"""
    FubiniStudy <: Metric
    rmse(x, y)

Fubini Study distance in complex projective plane PC^2
"""
struct FubiniStudy <: DD.Metric end
evaluate(::FubiniStudy, a, b) = _FubiniStudy(a, b)
(::FubiniStudy)(a, b) = _FubiniStudy(a, b)

get_distances(normal_nodes) = DD.pairwise(DIST[], normal_nodes)

######################
### Pairwise distances
######################
# lowertriangular_bool(I...) = LA.tril!(trues(I), -1)
# upperdiagonal(N::Int) = CC.combinations(1:N, 2)

# function pairwise_FS(nodes::AbstractArray, N::Int = length(nodes))
#     out = zeros(Float64, N, N)
#     f(I; fnodes = nodes) = _FubiniStudy(fnodes[I[1]], fnodes[I[2]])
#     I = lowertriangular_bool(N, N)
#     J = findall(I)
#     out[I] .= f.(J)
#     return LA.symmetric(out, :L)
# end

######################
### Pairwise distances
######################
# get_distance(dists, I...; dict = lt_indices(dists)) = dists[dict[I]]

# _eltype(A::T) where {T<:AbstractArray} = _eltype(first(A))
# _eltype(a) = typeof(a)

# sym_to_vec(M) = M[_lt_indices(size(M))]

# function vec_to_sym_size(N::Int)
#     n = fld(isqrt(8*N + 1) + 1, 2)
#     n * (n - 1) == 2 * N || @warn "vec_to_sym: length of vector is not triangular."
#     return n
# end
# vec_to_sym_size(v) = vec_to_sym_size(length(v))

# # vec_to_sym_mat(::Type{T}, n::Int) where T = zeros(T, n, n)
# vec_to_sym_mat(v) = (n = vec_to_sym_size(v); zeros(_eltype(v), n, n))

# function vec_to_sym(v)
#     sym = vec_to_sym_mat(v)
#     sym[_lt_indices(size(sym)...)] .= v
#     return LA.symmetric(sym, :L)
# end

# struct TriangularIndices{N} end
# function Base.iterate(ti::TriangularIndices{N}) where N
#     return N:N, (0,0,1)
# end
# function Base.iterate(ti::TriangularIndices{N}, I) where N
#     k, l, i = I
#     k += i; i += 1; l += i
#     !(l < N) && return nothing
#     out = (N-l):(N-k)
#     return out, (k, l, i)
# end

# function triangularindices_cols(N::Int)
#     out = UnitRange{Int64}[]
#     k, l, i = 0, 0, 1
#     while l < N
#         push!(out, (N-l):(N-k))
#         k += i; i += 1; l += i
#     end
#     return reverse(out)
# end
# triangularindices_cols(v) = triangularindices_cols(length(v))

# function triangularindices(x)
#     cols = triangularindices_cols(x)
#     k = 0
#     out = Vector{Tuple{Int, Int}}[]
#     for (l, I) in enumerate(cols)
#         J, k = I.-(k-l), last(I)
#         push!(out, tuple.(J, [l]))
#     end
#     return out
# end

# function lt_mat(v)
#     indices = triangularindices(v)
#     n = length(first(indices)) + 1
#     out = zeros(eltype(v), n, n)
#     for (i, key) in enumerate(reduce(vcat, indices))
#         out[key...] = v[i]
#     end
#     return out
# end
