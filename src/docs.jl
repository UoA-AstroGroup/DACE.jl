# add some documentation
@doc """
DA object representing a single polynomial.

The default constructor (with no arguments) creates an empty DA object
representing the constant zero function.

""" DA

@doc """
    DA(c::Real)

Create a DA object with the constant part equal to `c`.
""" DA(c::Real)

@doc """
    DA(i::Integer, c::Real)

Create a DA object as `c` times the independent variable number `i`.
Index zero creates the constant `c`.
""" DA(i::Integer, c::Real)

@doc """
    init(order, variables)

Initialize the global DACE engine. Both arguments must be positive. Existing DA
values become invalid, even if the dimensions are unchanged. Configure the
engine before launching concurrent calculations.
""" init

@doc """
    compile(p::DA)
    compile(polynomials::AbstractVector{<:DA})

Build an evaluation tree once for repeated numeric evaluation or DA composition.
`compiledDA` is the equivalent constructor. Inspect its output dimension, maximum
variable index, order and term count with `getDim`, `getVars`, `getOrd`, `getTerms`.
""" compile

@doc """
    evaluate(p, coordinates)
    evaluate(compiled, coordinates, result)

Evaluate a scalar polynomial, vector of polynomials, or compiled map. Real
coordinates produce numbers; DA coordinates compose polynomials. Missing
coordinates are zero and extra coordinates are ignored. Scalar polynomials
return scalars; maps return vectors. The three-argument form writes into `result`.

This replaces polynomial `DACE.eval`; Julia's module evaluator is unchanged.
""" evaluate

@doc """
    evaluate!(result::Vector{Float64}, compiled, coordinates, work::Vector{Float64})

Evaluate a compiled map numerically without allocating when all coordinates are
supplied. `result` must have `getDim(compiled)` entries and `work` must have at
least `getOrd(compiled)+1` entries. These buffers and the inputs must not alias.
Use `evaluate` for polynomial composition.
""" evaluate!

@doc """
    getCoefficient(p, exponents)
    getCoefficient(m::Monomial)

Read a Taylor coefficient (without factorial scaling). Missing exponents are
zero, extra exponents are ignored, and orders above the initialized maximum
return zero. Use `setCoefficient!` to change a coefficient in place.
""" getCoefficient

@doc """
    norm(p::DA, type=0)

Coefficient norm: type zero is the maximum absolute coefficient; type one is
their absolute sum; higher integer types use the corresponding power norm.
Unlike DA comparisons, this includes nonconstant coefficients.
""" norm

@doc """
    getRawMoments(mgf, order)

Return `(multi_indices, moments)` through the requested total order. `mgf` must
be a moment-generating function about zero. Taylor coefficients are multiplied
by the factorials of their exponents to recover the raw moments.
""" getRawMoments

@doc """
    getCentralMoments(mgf, order)

Return `(multi_indices, moments)` after centering the moment-generating function
at its mean. The input is assumed to have constant coefficient one.
""" getCentralMoments
