# Working with a fix partition.
# We store it in a callable struct 'Fitness'

function _debug_mpp(curve, multiplepoints, partition)
    @debug begin
        arenodes = check_multiplepoints(curve, multiplepoints)
        "Multiplepoints checked:", arenodes
    end
    @debug "Distances..."
    @debug begin
        M = get_distances(curve, multiplepoints)
        "Perimeter: ", total_perimeter(partition, M)
    end
end

struct Fitness
    curve::Vector{Float64}
    multiplepoints::Vector{Vector{ComplexF64}}
    partition::Vector{Vector{Int}}
    # F::T
    # Fitness(partition::Vector{Vector{Int}}) = new(partition)
end

function Fitness(curve::AbstractVector{<:Real} = randcurve(); kwargs...)
    @debug "Verbose debugging information.  Invisible by default"
    @debug "Computing fit function for curve."
    multiplepoints = get_multiplepoints(curve; kwargs...)
    @debug "Partition..."
    partition = get_partition(curve, multiplepoints)
    _debug_mpp(curve, multiplepoints, partition)
    return Fitness(curve, multiplepoints, partition)
end

function (f::Fitness)(curve_end::AbstractVector{<:Real}; kwargs...)
    # multiplepoints = track_multiplepoints(f.multiplepoints, f.F, curve; kwargs...)
    # multiplepoints = get_multiplepoints(curve; kwargs...)
    @debug "Verbose debugging information.  Invisible by default"
    @debug "Computing curve's fitness:"
    @debug "Traking multiplepoints..."
    @unpack curve, partition, multiplepoints = f
    multiplepoints_end = track_multiplepoints_flat(curve, curve_end, multiplepoints)
    _debug_mpp(curve_end, multiplepoints_end, partition)
    return total_perimeter(curve_end, multiplepoints_end, partition)
end

function Base.print(io::IO, f::Fitness, kwargs...)
    _print_partition(io, f.partition)
    !(isempty(kwargs)) && print(io, kwargs)
end
Base.show(io::IO, ::MIME"text/plain", f::Fitness) = Base.print(io, f)

struct FastFitness{T}
    multiplepoints::Vector{Vector{ComplexF64}}
    partition::Vector{Vector{Int}}
    param_homotopy::T
end

FastFitness(_fit::Fitness, param_homotopy) =
    FastFitness(_fit.multiplepoints, _fit.partition, param_homotopy)
function FastFitness(curve::AbstractVector{<:Real} = randcurve(); kwargs...)
    _fit = Fitness(curve; kwargs...) # return Fitness(curve, multiplepoints, partition)
    @debug "Param param system..."
    param_homotopy =
        param_param_system(complexfy(curve), PARAMS_END[], PARAM[]; gamma = randn())
    return FastFitness(_fit, param_homotopy)
end

function (ffit::FastFitness)(curve_end::AbstractVector{<:Real}; kwargs...)
    @debug "Verbose debugging information.  Invisible by default"
    @debug "Computing curve's fitness:"
    @debug "Traking multiplepoints..."
    @unpack multiplepoints, partition, param_homotopy = ffit
    multiplepoints_end =
        param_track_multiplepoints_flat(param_homotopy, curve_end, multiplepoints)
    _debug_mpp(curve_end, multiplepoints_end, partition)
    return total_perimeter(curve_end, multiplepoints_end, partition)
end

function Base.print(io::IO, ffit::FastFitness; kwargs...)
    _print_partition(io, ffit.partition)
    !(isempty(kwargs)) && print(io, kwargs)
end
Base.show(io::IO, ::MIME"text/plain", ffit::FastFitness) = Base.print(io, ffit)
