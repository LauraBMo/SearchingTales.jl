## Find the multiple points of a planar curve given by its parametritation (which is an array of polynomials).
## That is, given a parametritation of a curve as an array of polynomials,
## we get the multiple points (or the singular points) of the curve.

@doc raw"""
    diagonal_coeffs((P, Q))

Returns the coefficitens to the equations ... in a matrix.

Given variables x,y and polynomials w/o common roots:
P(x) = a_0 + a_1*x ... + a_n*x^n ≅ [a_0,a_1,..,a_n] = v
Q(x) = b_0 + b_1*x + ...+ b_n*x^n ≅ [b_0,b_1,..,b_n] = u
Then, entries of M are the coefficients of the polynomial:
P(x)Q(y) - P(y)Q(x) = sum(M[i,j]*x^(i-1)*y^(j-1), i, j)
That is, M == 0 iff [P(x) : Q(x)] = [P(y) : Q(y)]
or, ([P(x) : Q(x)], [P(y) : Q(y)]) ∈ diagonal ⊂ P1xP1

# Examples
```jldoctest
julia> diagonal_coeffs(curve_to_polys(ones(66))[1:2])
insert result of diagonal_coeffs(curve_to_polys(ones(66))[1:2])
```
"""
function diagonal_coeffs((v, u))
    M = v * LA.transpose(u)
    # Recover `vec(M)` with `isqrt ∘ length`
    return M - LA.transpose(M)
end

# function curve_to_multiplepoints_eq(curve)
#     # println("Polys: ", length(polynomials), ", polynomials' length: ", length.(polynomials))
#     diagonal_coeffs.(CC.combinations(curve_to_polys(curve), 2))
# end

# Build a generic dense system of polynomials suited for the polynomials
# Convert Julia Cartesian indices to polynomial support/exponents arrays.
# index_to_support(I) = (collect ∘ Tuple)(I) .- 1
# SupportIndices(axes::Tuple) = index_to_support.(CartesianIndices(axes))

@doc raw"""
    SupportIndices(A::AbstractArray)

Returns array of size `size(A)` whose `I` component is `I .-1`.

# Examples
```jldoctest
julia> SupportIndices((3,3))
3×3 Matrix{Vector{Int64}}:
 [0, 0]  [0, 1]  [0, 2]
 [1, 0]  [1, 1]  [1, 2]
 [2, 0]  [2, 1]  [2, 2]
```
"""
function SupportIndices end
SupportIndices(axes::Tuple) = collect.(Iterators.product([0:(i-1) for i in axes]...))
SupportIndices(A::AbstractArray) = SupportIndices(size(A))

function diagonal_polys(polys)
    coeffs = diagonal_coeffs(polys)
    exponents = SupportIndices(coeffs)
    out = vec(collect(zip(exponents, coeffs)))
    # sort!(out, lt=HC.ModelKit.td_order, by=first)
    return hstack(first.(out)), last.(out)
end

function diagonal_system(curve; kwargs...)
    @debug "Final curve:" length(curve_to_polys(curve))
    @debug "Combinations:" length(CC.combinations(curve_to_polys(curve), 2))
    supps_coeffs = diagonal_polys.(CC.combinations(curve_to_polys(curve), 2))
    supps, coeffs = first.(supps_coeffs), last.(supps_coeffs)
    # @debug "Coeffs: ", typeof(coeffs), ", lengths: ", length.(coeffs)
    # @debug "Supports: ", typeof(supps), ", matrix sizes: ", size.(supps)
    # @debug "Coefficients: ", coeffs
    return HC.System(supps, coeffs; variables=VARS[], kwargs...)
end

get_multiplepoints(curve, F = diagonal_system(curve); kwargs...) =
    sort_byreal!(_solve_onlynonsingular(F; kwargs...))


## Get multiple-points via Monodromy




## Track multiple-points of 'curve_init' to the ones of 'curve_end'
## Build flat homotopy H(x,param)
function param_system(curve_init, curve_end, param; gamma)
    @debug "Param curve: " length(curve_init) length(curve_end)
    param_curve = gamma .* (param .* curve_init) + (1 - param) .* curve_end
    return diagonal_system(param_curve; parameters=[param])
end

function track_multiplepoints_flat(curve_init, curve_end, multiplepoints;
                                   gamma=randn(),
                                   kwargs...)
    @debug "Verbose debugging information.  Invisible by default"
    @debug "Tracking multiplepoints (flat homotopy)"
    @debug "Setting homotopy..."
    # param = HC.Variable(gensym(:t))
    param = PARAM[]
    vars = VARS[]
    Gt = param_system(curve_init, curve_end, param; gamma=1.0)
    homotopy = HC.ParameterHomotopy(HC.fixed(Gt; compile=false), [1.0], [0.0])
    # homotopy = HC.ParameterHomotopy(Gt, [1.0], [0.0])
    # homotopy = HC.Homotopy(HC.expressions(Gt), vars, param)
    @debug "Solving..."
    result = _solve(homotopy, multiplepoints; kwargs...)
    # result = _solve(homotopy, multiplepoints;
    #     # seed=0x75a6a462,
    #     kwargs...)
    # println(result)
    # println(HC.solutions(result;
    #     only_real=false,
    #     real_tol=1e-6,
    #     only_nonsingular=false,
    #     only_singular=false,
    #     only_finite=false,
    #     multiple_results=false,
    # ))
    # out = HC.solutions(result;
    #     only_real=false,
    #     # real_tol=1e-6,
    #     only_nonsingular=true,
    #     only_singular=false,
    #     only_finite=true,
    #     multiple_results=false,
    # )
    out = HC.solutions(result)
    @debug "Solved! Solutions: $(length(out))"
    # return sort_byreal!(out)
    return out
end

# function track_multiplepoints(multiplepoints, F, curve_new; kwargs...)
#     G = diagonal_system(curve_new)
#     homotopy = HC.StraightLineHomotopy(G, F; gamma=randn())
#     result = _solve(homotopy, multiplepoints;
#                     # seed=0x75a6a462,
#                     kwargs...)
#     println(result)
#     # println(HC.solutions(result;
#     #     only_real=false,
#     #     real_tol=1e-6,
#     #     only_nonsingular=false,
#     #     only_singular=false,
#     #     only_finite=false,
#     #     multiple_results=false,
#     # ))
#     return HC.solutions(result;
#         only_real=false,
#         real_tol=1e-6,
#         only_nonsingular=false,
#         only_singular=false,
#         only_finite=true,
#         multiple_results=false,
#     )
# end

eval_nodes(curve, multiplepoints) =
    LA.normalize.(evalcurve.([curve], first.(multiplepoints)))

function get_nodes(curve, F=diagonal_system(curve); kwargs...)
    multiplepoints = get_multiplepoints(curve, F; kwargs...)
    return eval_nodes(curve, multiplepoints)
end

function check_multiplepoint(curve, coordinates; kwargs...)
    notindiag = !(isapprox(first(coordinates), last(coordinates)))
    pxs, pas = evalcurve.([curve], coordinates)
    # Are all the maximal minors px*qa-pa*qx zero?
    isoverlap = equal_projective(pxs, pas; kwargs...)
    return notindiag && isoverlap
end

# Check multiple nodes (namely, the 36 of a random planar curve of degree 10)
check_multiplepoints(curve, multiplepoints; kwargs...) =
    all(mp -> check_multiplepoint(curve, mp; kwargs...), multiplepoints)
