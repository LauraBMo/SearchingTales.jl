
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

struct Fitness
    curve::Vector{Float64}
    partition::Vector{Vector{Int}}
    multiplepoints::Vector{Vector{ComplexF64}}
    # F::T
    # Fitness(partition::Vector{Vector{Int}}) = new(partition)
end

function Fitness(curve::AbstractVector{<:Real} = randcurve(); kwargs...)
    F = diagonal_system(curve)
    @debug "Verbose debugging information.  Invisible by default"
    @debug "Computing fit function for curve."
    @debug "Getting multiple points..."
    multiplepoints = get_multiplepoints(curve, F; kwargs...)
    @debug "Eval nodes..."
    nodes = eval_nodes(curve, multiplepoints)
    @debug begin
        arenodes = check_multiplepoints(curve, multiplepoints)
        "Nodes checked: $(arenodes)"
    end
    @debug "Distances..."
    M = get_distances(nodes)
    return Fitness(curve, get_partition(M), multiplepoints)
end

function (f::Fitness)(curve::AbstractVector{<:Real}; kwargs...)
    # multiplepoints = track_multiplepoints(f.multiplepoints, f.F, curve; kwargs...)
    # multiplepoints = get_multiplepoints(curve; kwargs...)
    @debug "Verbose debugging information.  Invisible by default"
    @debug "Computing curve's fitness:"
    @debug "Traking multiplepoints..."
    multiplepoints = track_multiplepoints_flat(f.curve, curve, f.multiplepoints)
    @debug "Eval to get nodes..."
    nodes = eval_nodes(curve, multiplepoints)
    @debug begin
        arenodes = check_multiplepoints(curve, multiplepoints)
        "Nodes checked: $(arenodes)"
    end
    # println(_dehomo.(nodes))
    @debug "Distances..."
    M = get_distances(nodes)
    # return total_perimeter(f.partition, M), nodes
    @debug "Perimeter..."
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

function param_param_system(curve_init, curve_params, param; gamma)
    pp_curve = gamma .* (param .* curve_init) + (1 - param) .* curve_params
    @debug "Final curve:" length(splitDIM(pp_curve))
    @debug "Combinations:" length(CC.combinations(splitDIM(pp_curve), 2))
    supps_coeffs = diagonal_polys.(CC.combinations(splitDIM(pp_curve), 2))
    supps, coeffs = first.(supps_coeffs), last.(supps_coeffs)
    exprs = HC.horner.(build_system(supps, coeffs, VARS[]), [curve_params])
    homotopy = HC.Homotopy(exprs, VARS[], param, curve_params)
    # @debug "Coeffs: ", typeof(coeffs), ", lengths: ", length.(coeffs)
    # @debug "Supports: ", typeof(supps), ", matrix sizes: ", size.(supps)
    # @debug "Coefficients: ", coeffs
    # HC.System(supps, coeffs; variables=VARS[], kwargs...)
    # _parameters = curve_params; pushfirst!(_parameters, param)
    # return diagonal_system(param_param_curve; parameters=_parameters)
    return homotopy
end

struct FastFitness{T}
    partition::Vector{Vector{Int}}
    multiplepoints::Vector{Vector{ComplexF64}}
    param_homotopy::T
end

function FastFitness(curve::AbstractVector{<:Real} = randcurve(); kwargs...)
    # F = diagonal_system(curve)
    @debug "Verbose debugging information.  Invisible by default"
    @debug "Computing fit function for curve."
    @debug "Getting multiple points..."
    multiplepoints = get_multiplepoints(curve; kwargs...)
    @debug begin
        arenodes = check_multiplepoints(curve, multiplepoints)
        "Nodes checked: $(arenodes)"
    end
    @debug "Partition..."
    partition = get_partition(curve, multiplepoints)
    # tracker(newcurve) = track_multiplepoints_flat(curve, newcurve, multiplepoints)
    homotopy = param_param_system(complexfy(curve), PARAMS_END[], PARAM[]; gamma = randn())
    return FastFitness(partition, multiplepoints, homotopy)
end

function (ffit::FastFitness)(newcurve::AbstractVector{<:Real}; kwargs...)
    # multiplepoints = track_multiplepoints(f.multiplepoints, f.F, curve; kwargs...)
    # multiplepoints = get_multiplepoints(curve; kwargs...)
    @debug "Verbose debugging information.  Invisible by default"
    @debug "Computing curve's fitness:"
    @debug "Traking multiplepoints..."
    homo = HC.fix_parameters(ffit.param_homotopy, newcurve; compile = :all)
    ## :all, :none, :mixed
    newmultiplepoints = HC.solutions(_solve(homo, ffit.multiplepoints; kwargs...))
    @debug begin
        arenodes = check_multiplepoints(newcurve, newmultiplepoints)
        "Nodes checked: $(arenodes)"
    end
    # println(_dehomo.(nodes))
    @debug "Distances..."
    newnodes = eval_nodes(newcurve, newmultiplepoints)
    M = get_distances(newnodes)
    # return total_perimeter(f.partition, M), nodes
    @debug "Perimeter..."
    return total_perimeter(ffit.partition, M)
end

function Base.print(io::IO, ffit::FastFitness)
    print(io, "┌Partition:\n")
    l = 3; i = 1
    while i < l+1
        print(io, ffit.partition[i:l:end], "\n")
        i += 1
    end
    # !(isempty(kwargs)) && print(io, kwargs)
end

Base.show(io::IO, ::MIME"text/plain", ffit::FastFitness) = Base.print(io, ffit)
