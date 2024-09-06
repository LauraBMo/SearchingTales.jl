
# _nonsingularpaths(results) = filter!(HC.is_nonsingular, HC.path_results(results))
# unique_multiplepoints(points; atol) = uniquetol(p -> sort_byreal(p...), points; atol = atol)

# HC_multiplepoints(solutions; atol = 0) = unique_multiplepoints(solutions; atol > 0 ? atol : 1e2*accuracy)

# function fitness(curve::AbstractVector{<:Real};
#                  nodes = get_nodes(curve),
#                  M = get_distances(nodes),
#                  partition = get_partition(M),
#                  kwargs...)
#     # M = DD.pairwise(DIST[], reduce(hcat, nodes), dims = 2)
#     # println("Fit: Nodes, M, T computed!", nodes, " ", M, " ", partition)
#     # println("Fit: Nodes, M, T computed!", " ", M, " ")
#     # multiplepoints, partition = newcurve_params(curve; kwargs...)
#     # nodes, partition = get_info(DICT[], curve)
#     return total_perimeter(partition, M)
# end

# fitness(curve::AbstractVector{<:Real}, partition; kwargs...) = fitness(curve;
#                                                                        partition = partition,
#                                                                        kwargs...)

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

# Working with a fix partition.
# We store it in a callable struct 'Fitness'

struct Fitness{T}
    curve::Vector{Float64}
    partition::Vector{Vector{Int}}
    multiplepoints::Vector{Vector{ComplexF64}}
    F::T
    # Fitness(partition::Vector{Vector{Int}}) = new(partition)
end

function Fitness(curve::AbstractVector{<:Real}; kwargs...)
    F = diagonal_system(curve)
    multiplepoints = get_multiplepoints(curve, F; kwargs...)
    nodes = eval_nodes(curve, multiplepoints)
    M = get_distances(nodes)
    return Fitness(curve, get_partition(M), multiplepoints, F)
end

function (f::Fitness)(curve::AbstractVector{<:Real}; kwargs...)
    # multiplepoints = track_multiplepoints(f.multiplepoints, f.F, curve; kwargs...)
    # multiplepoints = get_multiplepoints(curve; kwargs...)
    multiplepoints = track_multiplepoints_flat(f.curve, curve, f.multiplepoints)
    nodes = eval_nodes(curve, multiplepoints)
    # println(_dehomo.(nodes))
    M = get_distances(nodes)
    # return total_perimeter(f.partition, M), nodes
    return total_perimeter(f.partition, M)
end

function Base.print(io::IO, f::Fitness, kwargs...)
    print(io, "┌Partition:\n")
    l = 3; i = 1
    while i < l+1
        print(io, f.partition[i:l:end], "\n")
        i += 1
    end
    !(isempty(kwargs)) && print(io, kwargs)
end

Base.show(io::IO, ::MIME"text/plain", f::Fitness) = Base.print(io, f)
