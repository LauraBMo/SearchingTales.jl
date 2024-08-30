
function hyper_mutate(rate = 0.6, gamma = 0.4)
    function mutate!(curve::AbstractVector; kwargs...)
        @debug "Mutate:"
        # old_curve = copy(curve)
        ## Mutation in-place
        for i in eachindex(curve)
            if rand() < rate
                curve[i] += randn(eltype(curve)) * gamma
            end
        end
        @debug "Mutated ok"
        # set_deformcurve!(DICT[], old_curve, curve; kwargs...)
        # newcurve_set!(DICT[], curve)
        @debug "Setting mutated ok"
        return curve
    end
    return mutate!
end

function hyper_crossover(a = () -> rand(), b = () -> 1. - a())
    function crossover(curve1::AbstractVector, curve2::AbstractVector; kwargs...)
        @debug "CrossOver:"
        child1 = a()*curve1 .+ b()*curve2
        child2 = b()*curve1 .+ a()*curve2
        @debug "Procreating ok"
        # set_deformcurve!(DICT[], curve2, child1; kwargs...)
        # set_deformcurve!(DICT[], curve1, child2; kwargs...)
        @debug "Setting child ok"
        return child1, child2
    end
    return crossover
end

##########################            ##########################
########################## Using Dict ##########################
##########################            ##########################

# function dict_fitness(curve::AbstractVector)
#     nodes, partition = get_info(curve)
#     fitness_perimeter(nodes, partition; distance = DIST[])
# end

# # const mp_type = Union{Nothing, Vector{ComplexF64}}
# const Tnode = Vector{ComplexF64}
# const Tdict = (Vector{ComplexF64}, Tuple{Vector{Tnode}, Vector{Int}})
# const DICT = Ref{Dict}(Dict{Tdict...}())

# function empty_dict!()
#     DICT[] = Dict{Tdict...}()
# end

# function get_info(key, dict = DICT[])
#     try
#         dict[key]
#     catch
#         @warn "Some key is missing."
#         newcurve_set!(dict, key)
#         dict[key]
#     end
# end

# function entry_curve!(curve::AbstractVector; kwargs...)
#     push!(DICT[], curve => curve_info(curve; kwargs...))
#     return curve
# end

# function curve_info(curve; kwargs...)
#     mp = get_multiplepoints(curve; kwargs...)
#     partition = get_partition(mp)
#     !(length(mp) == 36) && @error "Not 36 points!!!"
#     return mp, partition
# end

# # function newcurve_params(curve; kwargs...)
# #     mp = get_multiplepoints(curve; kwargs...)
# #     partition = get_partition(mp)
# #     !(length(mp) == 36) && @error "Not 36 points!!!"
# #     return mp, partition
# # end

# # function newcurve_set!(dict, curve::AbstractVector; kwargs...)
# #     push!(dict, curve => newcurve_params(curve; kwargs...))
# #     return curve
# # end

# # function set_curve_partition!(dict, curve::AbstractVector, multiplepoints)
# #     partition = get_partition(multiplepoints)
# #     delete!(dict, curve)
# #     push!(dict, curve => (multiplepoints, partition))
# #     return curve
# # end
# # set_curve_partition(curve::AbstractVector, multiplepoints) =
# #     set_curve_partition!(DICT[], curve, multiplepoints)

# # randcurve(::Type{T} = ComplexF64, l = 33; kwargs...) where T =
# #     newcurve_set!(DICT[], randn(T, l); atol = 0, kwargs...)

# function dict_reset!(dict, curve; kwargs...)
#     delete!(dict, curve)
#     push!(dict, key => get_multiplepoints(key))
# end
# dict_reset(curve::AbstractVector{<:Real}; kwargs...) = dict_reset!(DICT[], curve; kwargs...)

# function dict_get(dict, key; kwargs...)
#     try
#         dict[key]
#     catch
#         @warn "Some key is missing."
#         dict_reset!(dict, key; kwargs...)
#         dict[key]
#     end
# end
# dict_get(key; kwargs...) = dict_get(DICT[], key; kwargs...)

# function set_deformcurve!(dict, old_curve::AbstractVector, new_curve::AbstractVector; kwargs...)
#     prenodes = dict_get(dict, old_curve; kwargs...)
#     @debug "Get prenodes ok"
#     new_prenodes = track_nodes(prenodes, old_curve, new_curve)
#     @debug "Tracking ok"
#     I = findall(!isnothing, new_prenodes)
#     if length(I) < 36
#         @warn "Something $(36 - length(I)) went to nothing."
#         actual_mp = new_prenodes[I]
#         !(length(actual_mp) == 36) && @error "Not 36 points!!!"
#         # partition = get_partition(actual_mp)
#         push!(dict, new_curve => (actual_mp, partition))
#     else
#         push!(dict, new_curve => (new_prenodes, partition))
#     end
#     @debug "Push dict ok"
#     return old_curve
# end

# # function fitness(curve::AbstractVector)
# #     nodes, partition = get_info(DICT[], curve)
# #     fitness_perimeter(nodes, partition; distance = DIST[])
# # end

# function fitness(curve::AbstractVector; kwargs...)
#     multiplepoints, partition = newcurve_params(curve; kwargs...)
#     # nodes, partition = get_info(DICT[], curve)
#     nodes = _dehomo.(_evalpoly.([curve], first.(multiplepoints)))
#     # fitness_area(nodes, partition; distance = DIST[])
#     fitness_perimeter(nodes, partition; distance = DIST[])
# end
