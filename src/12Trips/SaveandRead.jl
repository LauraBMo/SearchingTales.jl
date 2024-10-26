__write(name, object) =
    open(name, "w") do io
        writedlm(io, object)
    end

## Save and read a curve and a fit function computing its current fitness.
_folder(name) = "12Trips/examples/" * name
function write(curve, fit::Fitness, name::AbstractString; folder = _folder(name))
    mkpath(folder)
    cd(() -> _pre_write(curve, fit), folder)
end
write(curve::AbstractVector{<:Real}, fit::Fitness, i::Int = 1; kwargs...) =
    write(curve, fit, "min$i"; kwargs...)

function _pre_write(curve, fit)
    multiplepoints = get_multiplepoints(fit, curve)
    ## Write info for `curve`
    __write("curve.txt", curve)
    __write("mp.txt", multiplepoints)
    perimeter = total_perimeter(curve, multiplepoints, fit.partition)
    @info "Perimeter of `curve`:" perimeter
    __write("fitness.txt", perimeter)
    mkpath("fit/")
    ## Write info for `fit`
    cd(() -> _pre_write_fit(fit), "fit/")
end

function _pre_write_fit(fit::Fitness)
    @unpack curve, partition, multiplepoints = fit
    __write("curve.txt", curve)
    __write("3part.txt", partition)
    __write("mp.txt", multiplepoints)
end
# Fitness{T}
#     curve::Vector{Float64}
#     multiplepoints::Vector{Vector{ComplexF64}}
#     partition::Vector{Vector{Int}}
#     param_homotopy::T

# Read
read(name::AbstractString; folder = _folder(name)) =
    cd(() -> _pre_read(), folder)
read(i::Int = 1; kwargs...) = read("min$i"; kwargs...)

function _pre_read()
    curve  = _pre_read_curve()
    # multiplepoints = _pre_read_vecvec("mp.txt", ComplexF64)
    old_fitness = first(readdlm("fitness.txt", Float64))
    fit = cd(() -> _pre_read_fit(), "fit/")
    new_fitness = fit(curve)
    # perimeter = total_perimeter(curve, multiplepoints, fit.partition)
    Delta = abs(old_fitness - new_fitness)
    @info "Fitness info:" old_fitness new_fitness Delta # perimeter
    return curve, fit
end

function _pre_read_fit()
    curve  = _pre_read_curve()
    partition = _pre_read_vecvec("3part.txt", Int)
    multiplepoints = _pre_read_vecvec("mp.txt", ComplexF64)
    @debug "Are multiple points:" check_multiplepoints(curve, multiplepoints)
    param_homotopy =
        param_param_system(curve; gamma = rand())
        # param_param_system(curve; gamma = 3.0)
    return Fitness(curve, multiplepoints, partition, param_homotopy)
end
# Fitness{T}
#     curve::Vector{Float64}
#     multiplepoints::Vector{Vector{ComplexF64}}
#     partition::Vector{Vector{Int}}
#     param_homotopy::T

function _pre_read_vecvec(name, ::Type{T}) where T
    out = readdlm(name, '\t', T, '\n')
    # Convert to convinient format
    return collect.(eachrow(out))
end
function _pre_read_curve()
    curve = readdlm("curve.txt", '\t', Float64, '\n')
    # Convert to convinient format
    return vec(curve)
end
