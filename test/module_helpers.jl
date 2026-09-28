# Ordinary modules use Core.EvalInto for eval on Julia 1.12+. DACE needs its own
# generic function so CxxWrap can add polynomial-evaluation methods to it.
@test DACE.eval(:(1 + 2)) == 3
@test DACE.eval(:(@__MODULE__)) === DACE

# Both include forms must evaluate in DACE, and the mapexpr form must actually
# apply its transformation. The fixture has no side effects on module globals.
fixture = joinpath(@__DIR__, "fixtures", "module_context.jl")
@test DACE.include(fixture) === DACE
module_name(expr) = expr isa LineNumberNode ? expr : :(nameof(@__MODULE__))
@test DACE.include(module_name, fixture) === :DACE

# Exercise the C++ overloads as well as the ordinary expression-eval helper.
DACE.init(3, 2)
x, y = DA(1, 1.0), DA(2, 1.0)
p = 1.0 + x + 2.0 * y + x * y
@test DACE.eval(p, AlgebraicVector([0.2, -0.3])) ≈ 0.54
q = DACE.eval(p, AlgebraicVector([y, x]))
@test DACE.eval(q, AlgebraicVector([0.2, -0.3])) ≈ 1.04
