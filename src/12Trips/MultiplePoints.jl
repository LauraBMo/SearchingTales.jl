## Find the multiple points of a planar curve given by its parametritation (which is an array of polynomials).
## That is, given a parametritation of a curve as an array of polynomials,
## we get the multiple points (or the singular points) of the curve.

@doc raw"""
    diagonal_coeffs((P, Q))

Returns the coefficitens to the equations ... in a matrix.

Given variables x,y and polynomials w/o common roots:
P(x) = a_0 + a_1*x ... + a_n*x^n ≅ [a_0,a_1,..,a_n]
Q(x) = b_0 + b_1*x + ...+ b_n*x^n ≅ [b_0,b_1,..,b_n]
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
    # out = Dict(exponents .=> coeffs_mat)
    sort!(out, lt=HC.ModelKit.td_order, by=first)
    return hstack(first.(out)), last.(out)
end

function diagonal_system(curve, parameters=nothing, variables=VARS[])
    # println(parameters, variables)
    supps_coeffs = diagonal_polys.(CC.combinations(curve_to_polys(curve), 2))
    supps, coeffs = first.(supps_coeffs), last.(supps_coeffs)
    # println("Coeffs: ", typeof(coeffs), ", lengths: ", length.(coeffs))
    # println("Supports: ", typeof(supps_coeffs))# ", matrix sizes: ", supps_coeffs)
    # println("Coefficients: ", coeffs)
    # println(parameters, variables)
    F =
        if isnothing(parameters)
            HC.System(supps, coeffs; variables=variables)
        else
            HC.System(supps, coeffs;
                      variables=variables,
                      parameters=parameters,
                      )
        end
    return F
end
# dense_support(axes::Tuple) = reduce(hcat, vec(SupportIndices(axes)))
# # dense_support(i::Int) = dense_support((i, i))
# dense_support(A::AbstractArray) = dense_support(size(A))


get_multiplepoints(curve, variables=VARS[], F=diagonal_system(curve, nothing, variables); kwargs...) =
   sort_byreal!(_solve_onlynonsingular(F; kwargs...))

function param_curve(curve_init, curve_end, param; gamma=randn())
    return gamma*(param.*curve_init) + (1-param)*curve_end
end

function track_multiplepoints_flat(curve_init, curve_end, multiplepoints; kwargs...)
    t = HC.Variable(gensym(:t))
    Ct = param_curve(curve_init, curve_end, t)
    Gt = diagonal_system(Ct, [t])
    homotopy = HC.ParameterHomotopy(HC.fixed(Gt; compile=false), [1.0], [0.0])
    result = _solve(homotopy, multiplepoints;
                    # seed=0x75a6a462,
                    kwargs...)
    # println(result)
    # println(HC.solutions(result;
    #     only_real=false,
    #     real_tol=1e-6,
    #     only_nonsingular=false,
    #     only_singular=false,
    #     only_finite=false,
    #     multiple_results=false,
    # ))
    out = HC.solutions(result;
        only_real=false,
        real_tol=1e-6,
        only_nonsingular=true,
        only_singular=false,
        only_finite=true,
        multiple_results=false,
    )
    return sort_byreal!(out)
end

function track_multiplepoints(multiplepoints, F, curve_new; kwargs...)
    G = diagonal_system(curve_new)
    homotopy = HC.StraightLineHomotopy(G, F; gamma=randn())
    result = _solve(homotopy, multiplepoints;
                    # seed=0x75a6a462,
                    kwargs...)
    println(result)
    # println(HC.solutions(result;
    #     only_real=false,
    #     real_tol=1e-6,
    #     only_nonsingular=false,
    #     only_singular=false,
    #     only_finite=false,
    #     multiple_results=false,
    # ))
    return HC.solutions(result;
        only_real=false,
        real_tol=1e-6,
        only_nonsingular=false,
        only_singular=false,
        only_finite=true,
        multiple_results=false,
    )
end

eval_nodes(curve, multiplepoints) =
    LA.normalize.(evalcurve.([curve], first.(multiplepoints)))

function get_nodes(curve, variables=VARS[], F=diagonal_system(curve, variables); kwargs...)
    multiplepoints = get_multiplepoints(curve, variables, F; kwargs...)
    return eval_nodes(curve, multiplepoints)
end

function check_multiplepoint(curve, coordinates; kwargs...)
    pxs, pas = evalcurve.([curve], coordinates)
    # Are all the maximal minors px*qa-pa*qx zero?
    return proj_eq(pxs, pas)
end

# Check multiple nodes (namely, the 36 of a random planar curve of degree 10)
check_multiplepoints(curve, nodes; kwargs...) =
    all(node -> check_multiplepoint(curve, node; kwargs...), nodes)
