## Some constants
## Number of VARiables, \P^1
const NVAR = 2
## Degree of the curve
const DEG = 10

# const CH_System = Ref{HC.System}()

## Dimension space of each parameter polynomial
const DIM = binomial(DEG+NVAR-1, NVAR-1) # 11
## Number of homogeneos coordinates in image
const N = 3
## Does not work:
# const X = Ref{Vector{HC.Variable}}([HC.Variable(:x), HC.Variable(:a)])
const VARS = Ref{Vector{HC.Variable}}()
const PARAM = Ref{HC.Variable}()
const PARAMS_END = Ref{Vector{HC.Variable}}()
const DIST = Ref{DD.Metric}()

function set_VARS(variables)
    VARS[] = variables
end

function set_PARAM(param)
    PARAM[] = param
end

function set_PARAMS_END(params)
    PARAMS_END[] = params
end

function set_metric(distance)
    DIST[] = distance
end

_multiexponents_affine(n, d) = Iterators.map(x -> (pop!(x); x), CC.multiexponents(n+1, d))

## Redefine for better performance.
function HC.ModelKit.monomials_exponents(n::Int, d::Int; affine::Bool = false)
    E =
        if affine
            _multiexponents_affine(n, d)
        else
            CC.multiexponents(n, d)
        end
    # HC.jl version, slow for large 'n' due to filtering
    # if affine
    #     pred = x -> sum(x) ≤ d
    # else
    #     pred = x -> sum(x) == d
    # end
    # E = map(Iterators.filter(pred, Iterators.product(Iterators.repeated(0:d, n)...))) do e
    #     collect(e)
    # end

    # println("n = ", n, " d = ", d)
    E = collect(E)
    sort!(E, lt = HC.ModelKit.td_order)
    E
end

function _solve_onlynonsingular(args...; kwargs...)
    # println("Solving multiple points...")
    # println("Args: ", args)
    nsolutions = 0
    top = 71 # The last solution must be already there.
    function areallfound(path)
        if HC.is_nonsingular(path)
            nsolutions += 1
        end
        return !(nsolutions < top)
    end
    @debug "Solving..."
    result = _solve(args...; stop_early_cb=areallfound, kwargs...)
    return HC.unique_points(HC.solutions(result); group_action = _flip)
end

function _solve_homotopy(homotopy, multiplepoints; kwargs...)
    @debug "Solving..."
    result = _solve(homotopy, multiplepoints; kwargs...)
    @debug result
    out = HC.solutions(result)
    @debug "Solved! nsolutions:", length(out)
    return out
end

# monomial(c, I, variables) = c * prod(variables .^ convert.(Int, I))
## Code from HomotopyContinuation.jl; model_kit/symbolic.jl::1128
# function build_system(support, coefficients, variables)
#     map(support, coefficients) do A, c
#         fi = HC.Expression(0)
#         for (k, I) in enumerate(eachcol(A))
#             HC.ModelKit.add!(fi, fi, monomial(c[k], I, variables))
#         end
#         fi
#     end
# end
# build_system(supp_coeffs, variables = VARS[]) =
#     build_system(first.(supp_coeffs), last.(supp_coeffs), variables)
# 6.936 ms (7294 allocations: 381.22 KiB)

__monomial(ci, ei) = ci * prod(VARS[] .^ convert.(Int, ei))
function __build_system(supp_coeffs)
    map(supp_coeffs) do (A, c)
        fi = HC.Expression(0)
        for (ci, ei) in zip(c, eachcol(A))
            HC.ModelKit.add!(fi, fi, __monomial(ci, ei))
        end
        fi
    end
end
# 7.005 ms (7294 allocations: 381.22 KiB)

# include("Utils.jl") ## Included in main
## Compute parameters for the multiple points of a parametric curve.
include("MultiplePoints.jl")
include("ThreePartition.jl")
include("Distances.jl")
set_metric(FubiniStudy())

include("Evolutionary.jl")

include("Fitness.jl")

## TODO
# include("TrackPoints.jl")
# include("Evolutionary.jl")


include("Recipes.jl")

include("SaveandRead.jl")
