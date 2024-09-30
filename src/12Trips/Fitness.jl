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

function Fitness(curve::AbstractVector{<:Real} = randcurve(), partition = nothing; kwargs...)
    @debug "Verbose debugging information.  Invisible by default"
    @debug "Computing fit function for curve."
    multiplepoints = get_multiplepoints(curve; kwargs...)
    if isnothing(partition)
        @debug "Partition..."
        partition = get_partition(curve, multiplepoints)
    end
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
function FastFitness(curve::AbstractVector{<:Real} = randcurve(), partition = nothing;
                     kwargs...)
    _fit = Fitness(curve, partition; kwargs...)
    # return Fitness(curve, multiplepoints, partition)
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

## Save and read a curve and its partition.
_folder(name) = "12Trips/examples/" * name
function _write(curve, partition, name; folder = _folder(name))
    mkpath(folder)
    cd((_) -> _pre_write(curve, partition), folder)
end
_write(curve, partition, i::Int = 1; kwargs...) =
    _write(curve, partition, "min$i"; kwargs...)

writeff(curve, ffit::FastFitness, name; kwargs...) =
    _write(curve, ffit.partition, name; kwargs...)
writeff(curve, ffit::FastFitness, i::Int = 1; kwargs...) =
    _write(curve, ffit, "min$i"; kwargs...)

function _pre_write(curve, partition)
    open("curve.txt", "w") do io
        writedlm(io, curve)
    end
    open("3part.txt", "w") do io
        writedlm(io, partition)
    end
end

# Read
_read(name::AbstractString; folder = _folder(name)) =
    cd(() -> _pre_read(), folder)
_read(i::Int = 1; kwargs...) = _read("min$i"; kwargs...)

function readff(name::AbstractString; kwargs...)
    curve, part = _read(name; kwargs...)
    # Create fast fitness funtion
    # from which to compute curve's fitness.
    return curve, FastFitness(randcurve(), part)
end
readff(i::Int = 1; kwargs...) = readff("min$i"; kwargs...)

function _pre_read()
    # Read data
    curve  = readdlm("curve.txt", '\t', Float64, '\n')
    part   = readdlm("3part.txt", '\t', Int, '\n')
    # Convert to convinient format
    curve = vec(curve)
    part = collect.(eachrow(part))
    return curve, part
end
