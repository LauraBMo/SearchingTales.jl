##
## Given a list of nodes `[1,2,...,n]` and its matrix of distances `M`,
## compute a 3-partition somehow minimizing the sum of the perimeters
## of such triangles.
##

# Compute the perimeter for a few random 3-partitions and take the minimum.
#
randpartition(N::Int) = collect.(Iterators.partition(Random.shuffle!(collect(1:N)), 3))

function rand_get_partition(curve::AbstractVector, multiplepoints, N = 20)
    M = get_distances(curve, multiplepoints)
    return rand_get_partition(M, 1:N)
end

function rand_get_partition(M::AbstractMatrix, I)
    dim = size(M, 1)
    @debug dim
    min_partition = randpartition(dim)
    init_perimeter = total_perimeter(min_partition, M)
    min_perimeter = init_perimeter
    for _ in I
        partition = randpartition(dim)
        perimeter = total_perimeter(partition, M)
        if perimeter < min_perimeter
            min_perimeter = perimeter
            min_partition = partition
            @debug min_perimeter
        end
    end
    @debug "gain", init_perimeter - min_perimeter
    return min_partition, min_perimeter
end

function rand_get_partition_multi(curve, multiplepoints, N = 20)
    M = get_distances(curve, multiplepoints)
    dim = size(M, 1)
    @debug dim
    I = 1:N
    chunks = Iterators.partition(I, length(I) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn rand_get_partition(M, chunk)
    end
    chunk_parts = fetch.(tasks)
    _, i = findmin(last, chunk_parts)
    return chunk_parts[i]
end

# Tring something smarter.
# Consider the complete graph over the 36 nodes.
# Sort edges by length, and start forming triangles
# from the shortest.
# So we get
# @btime PC.get_partition($C, $multiplepoints);
#   693.658 μs (22678 allocations: 1.54 MiB)
# To get a similar perimeter with the random methos above
# we need around  N = 1_000_000 with
# @btime PC.rand_get_partition($C, $multiplepoints, 1_000_000);
#   2.692 s (87000488 allocations: 5.51 GiB)

# edges:
# edge1 = i < j; edge2 = k < l
# pairs = [(i,k) (j,k))]
#         [(i,l) (j,l))]
allpairs(edges) = vec(collect(Iterators.product(edges...)))

findfirst_indiagonal(pairs) = findfirst(I -> reduce(==, I), pairs)

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
    get_partition(get_distances(curve, multiplepoints))

function _print_partition(io, partition)
    print(io, "┌ Partition:\n")
    l = 3; i = 1
    while i < l+1
        print(io, partition[i:l:end], "\n")
        i += 1
    end
    print(io, "└ End Partition\n")
end

function total_perimeter(triangles, M)
    loss = zero(eltype(M))
    for T in triangles
        # triangle = filter!(!isnothing, nodes[T])
        # if length(triangle) > 1 # it is a triangle or a line
        for (i, j) in CC.combinations(T, 2)
            loss += M[i, j]
        end
        # end # T is a point of less, so no loss increment.
    end
    return loss #, 1/loss
end

function total_perimeter(curve, multiplepoints, partition)
    M = get_distances(curve, multiplepoints)
    return total_perimeter(partition, M)
end
