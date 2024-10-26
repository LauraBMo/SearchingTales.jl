## Find the multiple points of a planar curve given by its parametritation (which is an array of polynomials).
## That is, given a parametritation of a curve as an array of polynomials,
## we get the multiple points (or the singular points) of the curve.

@doc raw"""
    diagonal_coeffs((P, Q))

Returns the coefficitens to the equations (see below) in a matrix.

Given variables `x,y` and polynomials w/o common roots:

`P(x) = a_0 + a_1*x ... + a_n*x^n ≅ [a_0,a_1,..,a_n] = v`

`Q(x) = b_0 + b_1*x + ...+ b_n*x^n ≅ [b_0,b_1,..,b_n] = u`

Then, entries of `M` are the coefficients of the polynomial:

`P(x)Q(y) - P(y)Q(x) = sum(M[i,j]*x^(i-1)*y^(j-1), i, j)`

That is, `M == 0` iff `[P(x) : Q(x)] = [P(y) : Q(y)]`,
or equivalently,\n

`([P(x) : Q(x)], [P(y) : Q(y)]) ∈ diagonal ⊂ P1 x P1`

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

@doc raw"""
    SupportIndices(A::AbstractArray)

Returns array of size `size(A)` whose `I`-th entry is `I.-1`.

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

"""
    diagonal_polys(polys)


"""
function diagonal_polys(polys)
    coeffs = diagonal_coeffs(polys)
    exponents = SupportIndices(coeffs)
    ## @btime fastfit(PC.randcurve())
    # #sort!(out, lt=HC.ModelKit.td_order, by=first)
    ## 7.756 ms (7288 allocations: 380.81 KiB)
    # sort!(out, lt=HC.ModelKit.td_order, by=first)
    ## 8.159 ms (7288 allocations: 380.81 KiB)
    ## Return pair:
    ## (matrix of exponents, vector of coeffitients)
    # out = zip(exponents, coeffs)
    # return hstack(first.(out)), last.(out)
    return hstack(exponents), coeffs
end

function _expressions(polys)
    combinations = CC.combinations(polys, 2)
    @debug "Final curve:" length(polys) length(combinations)
    ## Vector of pairs:
    ## (matrix of exponents, vector of coeffitients)
    supps_coeffs = diagonal_polys.(combinations)
    out = __build_system(supps_coeffs)                               ## [1]
    # out = HC.horner.(__build_system(supps_coeffs), [PARAMS_END[]]) ## [2]
    # out = HC.horner.(__build_system(supps_coeffs), [VARS[]])       ## [3]
    ## @btime fastfit(PC.randcurve())
    # [1] 8.622 ms (7288 allocations: 380.81 KiB)
    # [2] 8.035 ms (7291 allocations: 393.19 KiB)
    # [3] 8.487 ms (7291 allocations: 393.19 KiB)
    ## @btime fit(PC.randcurve())
    # [1] 6.731 s (8246247 allocations: 423.85 MiB)
    # [2] 6.778 s (8252456 allocations: 421.97 MiB)
    # [3] 6.474 s (7235093 allocations: 381.69 MiB)
    return out
end

function get_multiplepoints(curve; kwargs...)
    @debug "Multiple points..."
    exprs = _expressions(curve_to_polys(curve))
    F = HC.System(exprs; variables=VARS[])
    return _solve_onlynonsingular(F; kwargs...)
end

## Track multiple-points of 'curve_init' to the ones of 'curve_end'
## Build flat homotopy H(x,param)
##
function curve_homotopy(curve_init, curve_end = PARAMS_END[], t = PARAM[]; gamma)
    @debug "Param curve:" length(curve_init) length(curve_end)
    return gamma .* (t .* curve_init) + (1 - t) .* curve_end
end

function param_param_system(curve_init; gamma)
    params_param_curve = curve_homotopy(complexfy(curve_init); gamma)
    exprs = _expressions(splitDIM(params_param_curve))
    return HC.Homotopy(exprs, VARS[], PARAM[], PARAMS_END[])
end

function param_track_multiplepoints_flat(homotopy, curve_end, multiplepoints; kwargs...)
    fix_homo = HC.fix_parameters(homotopy, complexfy(curve_end); compile = :all)
    ## @btime fastfit(PC.randcurve())
    ## :all    #   8.685 ms (7288 allocations: 380.81 KiB)
    ## :mixed  # 126.428 ms (206692 allocations: 11.64 MiB)
    ## :none   # 136.295 ms (201992 allocations: 11.49 MiB)
    @debug "Traking multiplepoints..."
    return _solve_homotopy(fix_homo, multiplepoints; kwargs...)
end

function get_multiplepoints(fit, curve_end; kwargs...)
    @unpack multiplepoints, param_homotopy = fit
    multiplepoints_end =
        param_track_multiplepoints_flat(param_homotopy, curve_end, multiplepoints)
    @debug "Are multiple points:" check_multiplepoints(curve_end, multiplepoints_end)
    return multiplepoints_end
end

## Getting the actual singular points from the multiple points.
eval_nodes(curve, multiplepoints) =
    LA.normalize.(evalcurve.([curve], first.(multiplepoints)))

function get_nodes(curve; kwargs...)
    multiplepoints = get_multiplepoints(curve; kwargs...)
    return eval_nodes(curve, multiplepoints)
end

function check_multiplepoint(curve, coordinates; kwargs...)
    indiag = isapprox(first(coordinates), last(coordinates))
    indiag && @warn "Multiple point in diagonal!!!"
    pxs, pas = evalcurve.([curve], coordinates)
    # Are all the maximal minors px*qa-pa*qx zero?
    isoverlap = equal_projective(pxs, pas; kwargs...)
    return !(indiag) && isoverlap
end

# Check multiple nodes (namely, the 36 of a random planar curve of degree 10)
check_multiplepoints(curve, multiplepoints; kwargs...) =
    all(mp -> check_multiplepoint(curve, mp; kwargs...), multiplepoints)


# function param_system(curve_init, curve_end, t; gamma)
#     param_curve = curve_homotopy(curve_init, curve_end, t; gamma)
#     exprs = _expressions(curve_to_polys(param_curve))
#     return HC.System(exprs; variables=VARS[], parameters=[t])
# end

# function track_multiplepoints_flat(curve_init, curve_end, multiplepoints;
#                                    gamma=randn(),
#                                    kwargs...)
#     @debug "Verbose debugging information.  Invisible by default"
#     @debug "Tracking multiplepoints (flat homotopy)"
#     @debug "Setting homotopy..."
#     # t = HC.Variable(gensym(:t))
#     t = PARAM[]
#     Gt = param_system(curve_init, curve_end, t; gamma = gamma)
#     # homotopy = HC.ParameterHomotopy(HC.fixed(Gt; compile=true), [1.0], [0.0])
#     # homotopy = HC.ParameterHomotopy(HC.fixed(Gt; compile=false), [1.0], [0.0])
#     # homotopy = HC.ParameterHomotopy(Gt, [1.0], [0.0])
#     homotopy = HC.Homotopy(HC.expressions(Gt), VARS[], t)
#     return _solve_homotopy(homotopy, multiplepoints; kwargs...)
# end
