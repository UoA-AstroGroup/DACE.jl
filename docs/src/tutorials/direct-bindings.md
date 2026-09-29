# Direct bindings

DACE.jl calls DACE's C API with `ccall`. The native source is pinned to upstream
commit `368f6020342b0223ef1c05d50602bc54b7a8fe1f`. `deps/abi.c` only reports storage
layout and string-width information; it contains no polynomial algorithms.
CxxWrap, its generated allocated/reference types, the Julia interface fork, and
Eigen are absent from this implementation.

## API coverage

The scalar arithmetic, elementary and special functions, derivatives, integrals,
coefficient access, monomials, norms, bounds, truncation settings and factories
are retained. Compiled maps, numeric evaluation, DA composition, map inversion,
Jacobian/Hessian extraction, moments, and symmetric eigenpairs have Julia
implementations. `AlgebraicVector` retains componentwise arithmetic and
`AlgebraicMatrix` supports Julia linear algebra. Ordinary Julia arrays also work.

`DACE.evaluate` replaces polynomial `DACE.eval`, leaving Julia's module evaluator
alone. `evalScalar`, `compile` and `compiledDA` remain available. Missing evaluation
coordinates are zero; extra coordinates are ignored. Numeric compiled maps own
their coefficients and remain usable after engine reinitialization.

Comparisons retain the old constant-part semantics. Use `DACE.norm(p-q, 0)` to
compare all coefficients; `iszero(p)` also checks the entire polynomial. Powers
with a polynomial exponent now retain that exponent's derivatives. Julia's
`SpecialFunctions` functions work on DA values, including integer-order Bessel
functions. `multiplyMonomials` uses coefficient access to avoid an upstream routine
that fails to set its output length.

`eigh` computes real symmetric Taylor eigenpairs by solving each degree using
the constant bordered Jacobian. It supports distinct constant eigenvalues,
constant matrices, and diagonal polynomial matrices. Other repeated or numerically
indistinguishable eigenvalues are rejected: a unique analytic eigenbasis is not
generally defined there. The tests check fourth-order residuals and orthogonality,
as well as first derivatives against independent differentiation libraries.

## Ownership and concurrency

A `DA` owns its descriptor and coefficient buffer in Julia. `ccall` preserves both;
there are no native finalizers. Julia's collector accounts for the full buffer.
Constants and sparse arithmetic reserve only the space they can need;
`setCoefficient!` grows storage when necessary. `copy(p)` creates
independent storage; assignment shares the object, as with other mutable Julia
values. The Julia-backed array wrappers use the same reference semantics. Their
`copy` methods copy each scalar. `close(p)` optionally releases its buffer early.

`init` invalidates existing DA values and compiled-map composition. Such uses
throw an error instead of reading obsolete storage. Calls into the global C
engine are serialized with a reentrant lock, allowing tasks to migrate across
Julia threads. Numeric evaluation of compiled maps needs no native lock. Global
configuration (`init`, epsilon and truncation order) still affects all callers;
configure it before launching concurrent calculations.

## Performance

The DiffEqBase integration loads as a package extension when DiffEqBase is used;
ordinary DACE calculations do not load the SciML dependency stack.

Compile a map once when evaluating it repeatedly. Numeric `evaluate!` takes
caller-owned output and scratch buffers and allocates nothing when all needed
coordinates are supplied:

```julia
DACE.init(5, 2)
x, y = DACE.identity()
map = DACE.compile([sin(x)*exp(y), x+y*y])
args, result = [0.1, 0.2], zeros(2)
work = zeros(DACE.getOrd(map)+1)
DACE.evaluate!(result, map, args, work)
```

The input, output and scratch buffers must not alias. DA composition reuses
native multiplication and weighted-sum storage. Inversion solves one Taylor
degree at a time and fuses linear transformations to avoid quadratic numbers of
temporary polynomials. See the repository's `benchmark` folder for reproducible
comparisons with the previous CxxWrap implementation.
