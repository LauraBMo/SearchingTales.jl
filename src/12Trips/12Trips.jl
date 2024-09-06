## Some constants
## Number of VARiables, \P^1
const NVAR = 2
## Degree of the curve
const DEG = 10

# const CH_System = Ref{HC.System}()

## Dimension space of each parameter polynomial
const DIM = binomial(DEG+NVAR-1, NVAR-1)
## Number of homogeneos coordinates in image
const N = 3
## Does not work:
# const X = Ref{Vector{HC.Variable}}([HC.Variable(:x), HC.Variable(:a)])
const VARS = Ref{Vector{HC.Variable}}()

const DIST = Ref{DD.Metric}()

function set_VARS(variables)
    VARS[] = variables
end

function set_metric(distance)
    DIST[] = distance
end

multiexponents_affine(n, d) = Iterators.map(x -> (pop!(x); x), CC.multiexponents(n+1, d))

## Redefine for better performance.
function HC.ModelKit.monomials_exponents(n::Int, d::Int; affine::Bool = false)
    E =
        if affine
            multiexponents_affine(n, d)
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


using Colors
using RecipesBase
include("Recipes.jl")

using DelimitedFiles
include("SaveandRead.jl")
