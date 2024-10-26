module SearchingTales

using DocStringExtensions: SIGNATURES, TYPEDEF
using DelimitedFiles
using Colors
using RecipesBase
using Parameters

import HomotopyContinuation as HC
## My defaults for HC.jl (rightmost occurrence takes precedence in kwargs)
_solve(args...; kwargs...) = HC.solve(args...;
                                      start_system=:total_degree,
                                      show_progress=false,
                                      threading=false,
                                      kwargs...)

import LinearAlgebra as LA
import Combinatorics as CC
import Distances as DD
# using Evolutionary
using Random

include("Utils.jl")

include("12Trips/12Trips.jl")

end
