module DACE

using Libdl
using LinearAlgebra: LinearAlgebra, I, lu, SingularException
import SpecialFunctions

export DA, AlgebraicVector, AlgebraicMatrix, compiledDA, Monomial

include("engine.jl")
include("arithmetic.jl")
include("coefficients.jl")
include("arrays.jl")
include("evaluation.jl")
include("linear_algebra.jl")
include("statistics.jl")
include("docs.jl")

end
