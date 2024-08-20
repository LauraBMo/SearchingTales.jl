
## Save a curve and its info
get_curve(curve) = curve # Curve could be it own type.
function unpackall(curve)
    _curve = get_curve(curve) # convert to vector
    multiplepoints = get_multiplepoints(curve)
    nodes = eval_nodes(curve, multiplepoints)
    part = get_partition(curve, nodes)
    return _curve, multiplepoints, nodes, part
end
const EXTS = "_".*["curve", "multpp", "nodes", "part"].*".txt"
curve_iter(curve) = zip(EXTS, unpackall(curve))


function pre_save_PC(curve, name)
    for (naming, field) in curve_iter(curve)
        open(name*naming, "w") do io
            writedlm(io, field)
        end
    end
end
function save_PC(curve, name::AbstractString)
    mkpath(name)
    cd(() -> pre_save_PC(curve, name), name)
end

function pre_read_PC(name)
    names  = name.*EXTS
    # Read data
    curve  = readdlm(names[1], '\t', ComplexF64, '\n')
    multpp = readdlm(names[2], '\t', ComplexF64, '\n')
    nodes  = readdlm(names[3], '\t', ComplexF64, '\n')
    part   = readdlm(names[4], '\t', Int, '\n')
    # Convert to convinient format
    curve = vec(curve)
    multpp = Vector{ComplexF64}[collect(x) for x in eachrow(multpp)]
    nodes = Vector{ComplexF64}[collect(x) for x in eachrow(nodes)]
    part = collect.(eachrow(part))
    return curve, multpp, nodes, part
end
read_PC(name::AbstractString) = cd(() -> pre_read_PC(name), name)

save_PC(curve, i::Int = 1) = save_PC(curve, "minimizer$i")
read_PC(i::Int = 1) = read_PC("minimizer$i")
