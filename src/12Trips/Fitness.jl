# Working with a fix partition.
# We store it in a callable struct 'Fitness'

struct Fitness{T}
    curve::Vector{Float64}
    multiplepoints::Vector{Vector{ComplexF64}}
    partition::Vector{Vector{Int}}
    param_homotopy::T
end

function Fitness(curve::AbstractVector{<:Real} = randcurve(); kwargs...)
    @debug "Verbose debugging information. Invisible by default"
    multiplepoints = get_multiplepoints(curve; kwargs...)
    @debug "Are multiple points:" check_multiplepoints(curve, multiplepoints)
    partition = get_partition(curve, multiplepoints)
    param_homotopy =
        param_param_system(curve; gamma = rand())
        # param_param_system(curve; gamma = 1.0)
    return Fitness(curve, multiplepoints, partition, param_homotopy)
end

function (fit::Fitness)(curve_end::AbstractVector{<:Real}; kwargs...)
    @debug "Verbose debugging information.  Invisible by default"
    multiplepoints = get_multiplepoints(fit, curve_end; kwargs...)
    return total_perimeter(curve_end, multiplepoints, fit.partition)
end

## A few utils.jl
function Base.print(io::IO, f::Fitness, kwargs...)
    _print_partition(io, f.partition)
    !(isempty(kwargs)) && print(io, kwargs)
end
Base.show(io::IO, f::Fitness) = Base.print(io, f)
Base.show(io::IO, ::MIME"text/plain", f::Fitness) = Base.print(io, f)

## Same struct but ligther.
# struct FastFitness{T}
#     multiplepoints::Vector{Vector{ComplexF64}}
#     partition::Vector{Vector{Int}}
#     param_homotopy::T
# end

# FastFitness(fit::Fitness) =
#     FastFitness(fit.multiplepoints, fit.partition, fit.param_homotopy)
# FastFitness(curve::AbstractVector{<:Real} = randcurve(),
#             partition = nothing;
#             kwargs...) = FastFitness(Fitness(curve, partition; kwargs...))

# function (ffit::FastFitness)(curve_end::AbstractVector{<:Real}; kwargs...)
#     @debug "Verbose debugging information.  Invisible by default"
#     @debug "Computing curve's fitness:"
#     @debug "Traking multiplepoints..."
#     @unpack multiplepoints, partition, param_homotopy = ffit
#     multiplepoints_end =
#         param_track_multiplepoints_flat(param_homotopy, curve_end, multiplepoints)
#     # _debug_mpp(curve_end, multiplepoints_end, partition)
#     return total_perimeter(curve_end, multiplepoints_end, partition)
# end

# function Base.print(io::IO, ffit::FastFitness; kwargs...)
#     _print_partition(io, ffit.partition)
#     !(isempty(kwargs)) && print(io, kwargs)
# end
# Base.show(io::IO, ::MIME"text/plain", ffit::FastFitness) = Base.print(io, ffit)
