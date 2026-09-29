using DACE
using Test

include("utils.jl")

@testset verbose = true "DACE tests" begin
    @testset "Polynomial evaluation" begin
        DACE.init(3, 2)
        x, y = DA(1, 1.0), DA(2, 1.0)
        p = 1.0 + x + 2.0 * y + x * y
        @test DACE.eval(:(1 + 2)) == 3
        @test DACE.evaluate(p, AlgebraicVector([0.2, -0.3])) ≈ 0.54
        q = DACE.evaluate(p, AlgebraicVector([y, x]))
        @test DACE.evaluate(q, AlgebraicVector([0.2, -0.3])) ≈ 1.04
        for f in ([p, x + y], AlgebraicVector([p, x + y]), DACE.compile([p, x + y]))
            @test collect(DACE.evaluate(f, [0.2, -0.3])) ≈ [0.54, -0.1]
            g = DACE.evaluate(f, [y, x])
            @test collect(DACE.evaluate(g, [0.2, -0.3])) ≈ [1.04, -0.1]
        end
        result = AlgebraicVector(zeros(2))
        DACE.evaluate(DACE.compile([p, x + y]), AlgebraicVector([0.2, -0.3]), result)
        @test collect(result) ≈ [0.54, -0.1]
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
