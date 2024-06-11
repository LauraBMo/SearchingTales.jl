## Find the multiple points of a planar curve given by its parametritation (which is an array of polynomials).
## That is, given a parametritation of a curve as an array of polynomials,
## we get the multiple points (or the singular points) of the curve.


## Given variables x,y and polynomials w/o common roots:
## P(x) = a_0 + a_1*x ... + a_n*x^n ≅ [a_0,a_1,..,a_n]
## Q(x) = b_0 + b_1*x + ...+ b_n*x^n ≅ [b_0,b_1,..,b_n]
## Then, entries of M are the coefficients of the polynomial:
## P(x)Q(y) - P(y)Q(x) = sum(M[i,j]*x^(i-1)*y^(j-1), i, j)
## That is, M == 0 iff [P(x) : Q(x)] = [P(y) : Q(y)]
## or, ([P(x) : Q(x)], [P(y) : Q(y)]) ∈ diagonal ⊂ P1xP1
function diagonal_eq((P, Q))
    M = P * LA.transpose(Q)
    # Recover M with `isqrt ∘ length`
    return M - LA.transpose(M)
end
function curve_to_multiplepoints_eq(curve)
    # println("Polys: ", length(polynomials), ", polynomials' length: ", length.(polynomials))
    diagonal_eq.(CC.combinations(curve_to_polys(curve), 2))
end

# Build a generic dense system of polynomials suited for the polynomials
# Convert Julia Cartesian indices to polynomial support/exponents arrays.
index_to_support(I) = (collect ∘ Tuple)(I) .- 1
SupportIndices(axes) = index_to_support.(CartesianIndices(axes))
dense_support(axes::Tuple) = reduce(hcat, vec(SupportIndices(axes)))
# dense_support(i::Int) = dense_support((i, i))
dense_support(A::AbstractArray) = dense_support(size(A))


function get_multiplepoints(curve, variables=VARS[]; atol = 0, kwargs...)
    mp_coefficients = curve_to_multiplepoints_eq(curve)
    # println("Coeffs: ", typeof(mp_coefficients), ", its lengths: ", length.(mp_coefficients))
    # println("Supports: ", typeof(supp), ", matrix sizes: ", size.(supp))
    # println("Coefficients: ", mp_coefficients)
    supp = dense_support.(mp_coefficients)
    system = HC.System(supp,
        vec.(mp_coefficients);
        variables=variables,
    )
    return _solve_onlynonsingular(system; kwargs...)
end

get_nodes(curve, multiplepoints; atol = 0, kwargs...) =
   LA.normalize.(evalcurve.([curve], first.(multiplepoints)))

function get_nodes(curve; kwargs...)
    multiplepoints = get_multiplepoints(curve; kwargs...)
    return get_nodes(curve, multiplepoints)
end

function check_multiplepoint(curve, coordinates; kwargs...)
    pxs, pas = evalcurve.([curve], coordinates)
    # Are all the maximal minors px*qa-pa*qx zero?
    return proj_eq(pxs, pas)
end

# Check multiple nodes (namely, the 36 of a random planar curve of degree 10)
check_multiplepoints(curve, nodes; kwargs...) =
    all(node -> check_multiplepoint(curve, node; kwargs...), nodes)
