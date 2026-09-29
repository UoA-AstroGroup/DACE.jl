using DACE
using Test

include("utils.jl")

@testset verbose = true "DACE tests" begin
    @testset "Module and polynomial evaluation" begin
        @test DACE.eval(:(@__MODULE__)) === DACE
        DACE.init(3, 2)
        x, y = DA(1, 1.0), DA(2, 1.0)
        p = 1.0 + x + 2.0 * y + x * y
        @test DACE.eval(p, AlgebraicVector([0.2, -0.3])) ≈ 0.54
        q = DACE.eval(p, AlgebraicVector([y, x]))
        @test DACE.eval(q, AlgebraicVector([0.2, -0.3])) ≈ 1.04
    end

    @testset verbose = true "Tutorials" begin
        include("tutorial_tests.jl")
    end

    @testset verbose = true "Validation tests" begin
        include("validation_1.jl")
        include("validation_2.jl")
    end

    @testset verbose = true "Operators" begin
        include("comparison_operators.jl")
    end

    @testset verbose = true "Special Functions" begin
        include("special_functions.jl")
    end

    @testset verbose = true "Linear Algebra" begin
        include("linear_algebra.jl")
    end

    @testset verbose = true "Access & extraction" begin
        include("extraction.jl")
    end

    @testset verbose = true "Statistics" begin
        include("statistics.jl")
    end

    @testset verbose = true "Factories" begin
        include("factory.jl")
    end
end
