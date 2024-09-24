##
## Given a list of nodes `[1,2,...,n]` and its matrix of distances `M`,
## compute a 3-partition somehow minimizing the sum of the perimeters
## of such triangles.
##

findfirst_indiagonal(pairs) = findfirst(I -> reduce(==, I), pairs)

# edges:
# edge1 = i < j; edge2 = k < l
# pairs = [(i,k) (j,k))]
#         [(i,l) (j,l))]
allpairs(edges) = vec(collect(Iterators.product(edges...)))
function set_angle(i, pairs)
    common_node = first(pairs[i])
    sides = extrema(pairs[end-(i-1)])
    return [sides..., common_node]
end
function push_ifisangle!(angles, edges)
    pairs = allpairs(edges)
    i = findfirst_indiagonal(pairs)
    if !isnothing(i)
        push!(angles, set_angle(i, pairs))
    end
end

# 'map' runs over the 'zip' the iterators!
# angles and edges are sorted so...
istriangle(angle, edge) = all(map(==, angle, edge))

function findfirst_triangle(edges)
    angles = Vector{Int}[]
    for (ei, e) in enumerate(edges)
        i = findfirst(a -> istriangle(a, e), angles)
        if isnothing(i) # 'e' does _not_ form a triangle.
            for f in Iterators.take(edges, ei-1)
                push_ifisangle!(angles, (e, f))
            end
        else # 'e' _does_ form a triangle.
            return angles[i] # Angles and triangles are equal.
        end
    end
    return nothing
end

function get_partition!(edges)
    triangles = Vector{Int}[]
    while !isempty(edges)
        new_T = findfirst_triangle(edges)
        push!(triangles, new_T)
        # Remove all edges with some vertex in 'new_T'.
        filter!(edge -> isdisjoint(edge, new_T), edges)
    end
    return triangles
end

upperdiagonal(N::Int) = CC.combinations(1:N, 2)
function sorted_edges(M::AbstractMatrix)
    J = upperdiagonal(size(M, 1))
    dict = [Tuple(j) => M[j] for j in J]
    return first.(sort(dict, by = last))
end

get_partition(M::AbstractMatrix) = get_partition!(sorted_edges(M))
get_partition(curve, multiplepoints) =
    get_partition(get_distances(eval_nodes(curve, multiplepoints)))
