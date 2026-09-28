using DACE
using Test

include("utils.jl")

@testset verbose = true "DACE tests" begin
    @testset "Module helpers and polynomial evaluation" begin
        include("module_helpers.jl")
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
