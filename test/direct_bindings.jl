using LinearAlgebra, Random, SpecialFunctions

# DA equality intentionally compares constants. Check every coefficient instead.
polyerror(a::DA, b::Real) = DACE.norm(a-b, 0)
polyerror(a::AbstractArray, b::AbstractArray) = maximum(polyerror.(a, b); init=0.0)
@noinline temporary_buffer(p) = WeakRef(sin(p).buffer)

# A valid AbstractVector whose elements have no owner until materialized.
struct GeneratedPolynomials <: AbstractVector{DA}
    n::Int
end
Base.size(v::GeneratedPolynomials) = (v.n,)
function Base.getindex(v::GeneratedPolynomials, i::Int)
    checkbounds(v, i)
    GC.gc(false)
    DA(i,1.0)^2
end

@testset "Direct C bindings" begin
    @test DACE.eval(:(1+2)) == 3
    @testset "Constructors, coefficients and ownership" begin
        DACE.init(5, 3)
        DACE.setEps(0.0)
        x, y, z = DACE.identity()
        @test DACE.cons(DA(2)) == 2.0
        @test DACE.cons(DA(0, 2.5)) == 2.5
        @test DACE.cons(DA(true)) == 1.0
        @test DACE.cons(DA(1//2)) == 0.5
        p = 2 + 3x + 4x*y + 5z^3
        @test DACE.getMaxMonomials() == binomial(8, 3)
        @test DACE.getCoefficient(p, [1]) == 3
        @test DACE.getCoefficient(p, [1,1,0,9]) == 4
        @test DACE.getCoefficient(p, [6]) == 0
        @test DACE.linear(p) == [3,0,0]
        q = copy(p)
        DACE.setCoefficient!(q, [1,1], 7)
        @test DACE.getCoefficient(p, [1,1]) == 4
        @test DACE.getCoefficient(q, [1,1]) == 7
        ms = DACE.getMonomials(p)
        @test DACE.order.(ms) == [0,1,2,3]
        @test DACE.getCoefficient.(ms) == [2,3,4,5]
        @test DACE.getExponents(last(ms)) == [0,0,3]
        @test DACE.getCoefficient(DACE.getMonomial(p, 1)) == 2
        @test occursin("COEFFICIENT", DACE.toString(p))
        @test occursin("COEFFICIENT", sprint(show, p))
        @test !isempty(DACE.toString(DA()))
        @test !isempty(DACE.toString(first(ms)))
        @test_throws BoundsError DACE.getMonomial(p, 0)
        @test_throws ArgumentError DACE.setCoefficient!(p, [-1], 2)
        @test_throws ArgumentError DA(4, 1.0)
        @test_throws ArgumentError DACE.variable(0)
        close(q); close(q)
        @test_throws ArgumentError sin(q)
        @test isempty(q.buffer)
        # The descriptor and its entire backing buffer must stay alive across GC.
        for _ in 1:3
            GC.gc(true)
            @test polyerror(p, 2 + 3x + 4x*y + 5z^3) == 0
        end
        weakbuffer = temporary_buffer(p)
        GC.gc(true)
        @test weakbuffer.value === nothing
        @test Base.summarysize(p) >= sizeof(p.buffer)
        grown = DA()
        # Exercise buffer growth and pointer refresh, including across GC.
        indices = DACE.getMultiIndices(5,3)
        for (i,jj) in enumerate(indices)
            DACE.setCoefficient!(grown,jj,i)
            i % 7 == 0 && GC.gc(false)
        end
        @test DACE.getCoefficient.(Ref(grown),indices) == collect(1:length(indices))
        @test polyerror(copy(grown),grown) == 0
        compiled = DACE.compile(p)
        DACE.init(2, 1)
        @test_throws ArgumentError p+1
        @test_throws ArgumentError DACE.evaluate(compiled, [DA(1,1.0)])
        # Numeric evaluation owns its data and does not depend on engine state.
        @test only(DACE.evaluate(compiled, [0.1,0.2,0.3])) ≈ 2.515
        @test_throws ArgumentError DACE.init(0, 1)
        @test_throws ArgumentError DACE.init(typemax(Int), 2)
        @test DACE.cons(DA(7)) == 7 # Invalid init arguments leave the engine intact.
    end

    @testset "Scalar API and native error recovery" begin
        DACE.init(5, 2)
        DACE.setEps(0.0)
        x, y = DACE.identity()
        p = 2+x+0.2y
        for a in (2, 2.0, 2//1, π, ℯ)
            @test polyerror(a^(1+x), exp(log(a)*(1+x))) < 1e-14
        end
        @test polyerror(p^(1+y), exp(log(p)*(1+y))) < 1e-14
        @test polyerror(p^(1//2), sqrt(p)) < 1e-14
        @test polyerror(DACE.powi(p, 3), p^3) == 0
        @test polyerror(DACE.powd(p, 0.5), sqrt(p)) < 1e-14
        @test polyerror(DACE.isrt(p), inv(sqrt(p))) < 1e-14
        @test polyerror(DACE.icrt(p), inv(cbrt(p))) < 1e-14
        @test polyerror(log(3.0, p), log(p)/log(3.0)) < 1e-14
        @test polyerror(log(ℯ, p), log(p)) == 0
        @test polyerror(atan(x, p), atan(x/p)) < 1e-14
        @test polyerror(hypot(x,p), sqrt(x*x+p*p)) < 1e-14
        @test polyerror(mod(p+3, 3), p) == 0
        @test DACE.cons(round(DA(2.3)+x)) == 2
        @test DACE.cons(trunc(DA(-2.3)+x)) == -2
        @test iszero(sqrt(DA())) && iszero(cbrt(DA()))
        @test iszero(hypot(DA(), DA()))
        @test isless(1.0, p) && !isless(p, 1.0)
        @test isequal(2.0, p) && isequal(p, 2.0)
        @test (p == π) == (π == p) == false
        @test !iszero(x) && x == 0 # Legacy comparisons use the constant part.
        @test !iszero(DA(NaN)) && !iszero(DA(Inf))
        @test iszero(p-p)
        @test_throws DACE.DACEError log(-p)
        @test_throws DACE.DACEError inv(x)
        @test_throws DACE.DACEError sqrt(x)
        @test_throws ArgumentError p^typemin(Cint)
        @test polyerror(exp(log(p)), p) < 1e-13
        @test polyerror(DACE.deriv(p^3, [2,1]), 1.2) < 1e-13
        @test polyerror(DACE.deriv(DACE.integ(p, 1), 1), p) < 1e-14
        @test polyerror(DACE.deriv(DACE.integ(p, [1,1]), [1,1]), p) < 1e-14
        @test polyerror(DACE.divide(x^2*y, 1, 2), y) == 0
        @test polyerror(DACE.multiplyMonomials(2+3x, 4+5x+x*x), 8+15x) == 0
        @test DACE.norm(2+3x-4x*y, 0) == 4
        @test DACE.norm(2+3x-4x*y, 1) == 9
        @test DACE.orderNorm(2+3x-4x*y, 0, 1) == [2,3,4,0,0,0]
        @test length(DACE.estimNorm(exp(x), 0, 1, 8)) == 9
        @test_logs (:warn, r"estimate") DACE.estimNorm(x)
        @test iszero(DACE.multiplyMonomials(x, y))
        @test iszero(DACE.multiplyMonomials(x, DA()))
        # Sparse outputs reserve only as many terms as the operation can create.
        @test polyerror((2+x)*(3+y), 6+3x+2y+x*y) == 0
        @test polyerror((2+x)/DA(2), 1+0.5x) == 0
        @test polyerror(DA()-x, -x) == 0
        bounds = DACE.bound(2+3x)
        @test bounds.m_lb <= -1 && bounds.m_ub >= 5
        @test isfinite(DACE.random(0.2))
        @test DACE.getEpsMac() ≈ eps(Float64)
        @test DACE.setEps(1e-20) == 0
        @test DACE.getEps() == 1e-20
        @test DACE.setTO(3) == 5
        DACE.pushTO(2); DACE.pushTO(1)
        @test DACE.getTO() == 1
        DACE.popTO(); DACE.popTO()
        @test DACE.getTO() == 3
        @test_throws ArgumentError DACE.popTO()
        @test_throws ArgumentError DACE.setTO(6)
        DACE.setTO(5)
    end

    @testset "Special function derivatives" begin
        DACE.init(4, 1)
        DACE.setEps(0.0)
        x = 2.3 + DA(1,1.0)
        for f in (erf, erfc, gamma, loggamma, digamma)
            @test DACE.cons(f(x)) ≈ f(2.3) rtol=2e-13
        end
        for f in (besselj, bessely, besseli, besselk, besselix, besselkx), n in (0,1,2)
            @test DACE.cons(f(n,x)) ≈ f(n,2.3) rtol=2e-13
        end
        @test polyerror(DACE.trim(DACE.deriv(loggamma(x),1)-digamma(x),0,3), 0) < 1e-12
        @test polyerror(DACE.trim(DACE.deriv(besselj(1,x),1)-(besselj(0,x)-besselj(2,x))/2,0,3), 0) < 1e-12
        @test polyerror(DACE.PsiFunction(x,1), polygamma(1,x)) == 0
    end

    @testset "Compiled maps and inversion" begin
        DACE.init(5, 3)
        DACE.setEps(0.0)
        x,y,z = DACE.identity()
        f = [2+x+2y+x*z, z^3+y*y]
        compiled = compiledDA(f)
        @test (DACE.getDim(compiled), DACE.getVars(compiled), DACE.getOrd(compiled)) == (2,3,3)
        @test DACE.getTerms(compiled) > 1
        @test DACE.evaluate(compiled, [0.1,0.2,0.3]) ≈ [2.53,0.067]
        @test DACE.evaluate(compiled, [0.1]) ≈ [2.1,0.0]
        @test DACE.evaluate(compiled, [0.1,0.2,0.3,99]) ≈ [2.53,0.067]
        @test DACE.evalScalar(f, 0.1) ≈ [2.1,0.0]
        @test DACE.evalScalar(x, 0.3) == 0.3
        @test only(DACE.evaluate(DACE.compile(z), [1,2,3])) == 3
        @test DACE.evaluate(DACE.compile(GeneratedPolynomials(3)), [1,2,3]) == [1,4,9]
        @test DACE.evaluate(DA(7), Float64[]) == 7
        @test polyerror(DACE.evaluate(compiled, [x,y,z]), f) == 0
        @test polyerror(DACE.evaluate(compiled, Real[x,0,0]), [2+x,DA()]) == 0
        out, work = zeros(2), zeros(4)
        @test DACE.evaluate!(out,compiled,[0.1,0.2,0.3],work) ≈ [2.53,0.067]
        @test DACE.evaluate(compiled,[0.1,0.2,0.3],out) === out
        @test_throws DimensionMismatch DACE.evaluate!(zeros(1),compiled,[1.,2.,3.],work)
        @test_throws ArgumentError DACE.evaluate!(out,compiled,out,work)
        @test_throws ArgumentError DACE.evaluate!(out,compiled,[x,y,z],work)
        @test_throws ArgumentError DACE.compile(DA[])
        # A coupled map with a nontrivial Jacobian and a parameter coordinate.
        f = [2x+y+0.1x*y, x+3y+0.2x*x, z+0.1x*z]
        inverse = DACE.invert(f)
        @test polyerror(DACE.evaluate(f,inverse), [x,y,z]) < 1e-12
        @test polyerror(DACE.evaluate(inverse,f), [x,y,z]) < 1e-12
        partial = [x+0.2x*y+z*z, y+0.1x*z]
        inverse = DACE.invert(partial)
        @test polyerror(DACE.evaluate(partial,vcat(inverse,z)), [x,y]) < 1e-12
        @test_throws SingularException DACE.invert([x,x,z])
        @test DACE.getTO() == 5
        DACE.init(1,1)
        x = DA(1,1.0)
        @test polyerror(only(DACE.invert([2+3x])), (x-2)/3) < 1e-15
    end

    @testset "Arrays and higher-order symmetric eigenpairs" begin
        DACE.init(4,2)
        DACE.setEps(1e-24)
        x,y = DACE.identity()
        v = AlgebraicVector([1+x,2+y])
        @test polyerror(v*v, collect(v).^2) == 0
        @test polyerror(2/v, 2 ./collect(v)) == 0
        @test polyerror(sin(v), sin.(collect(v))) == 0
        @test_throws ArgumentError AlgebraicVector{DA}(-1)
        @test_throws ArgumentError AlgebraicMatrix{DA}(-1,2)
        filled = AlgebraicMatrix{DA}(2,2, x)
        DACE.setCoefficient!(filled[1,1], [0,0], 5)
        @test DACE.cons(filled[2,2]) == 0
        M = AlgebraicMatrix([2+x y; y 4-x])
        @test polyerror(DACE.det(M), (2+x)*(4-x)-y*y) < 1e-14
        @test polyerror(Matrix(M)*Matrix(inv(M)), DA.(Matrix{Float64}(I,2,2))) < 1e-13
        @test polyerror(Matrix(transpose(M)), Matrix(M)) == 0
        @test polyerror(DACE.frobenius(M), sqrt(sum(abs2,M))) == 0
        A = [2+x 0.3+y 0.2x; 0.3+y 4-y 0.1+x*y; 0.2x 0.1+x*y 7+x+y]
        values, vectors = DACE.eigh(AlgebraicMatrix(A))
        @test polyerror(A*vectors, vectors*Diagonal(values)) < 2e-12
        @test polyerror(transpose(vectors)*vectors, DA.(Matrix{Float64}(I,3,3))) < 2e-12
        @test DACE.getTO() == 4
        # Evaluate the Taylor eigenvalues against independent numeric LAPACK.
        for point in ([1e-3,2e-3], [-2e-3,1e-3])
            numeric = [DACE.evaluate(a,point) for a in A]
            @test DACE.evaluate(values,point) ≈ eigvals(Symmetric(numeric)) atol=2e-13
        end
        repeated, basis = DACE.eigh([1+x DA(); DA() 1+y])
        @test polyerror(repeated, [1+x,1+y]) == 0
        @test polyerror(basis, DA.(Matrix{Float64}(I,2,2))) == 0
        @test_throws ArgumentError DACE.eigh([1+x y; y 1-x])
        @test_throws ArgumentError DACE.eigh([x y; x y])
        @test size(DACE.hess_stack([x*x+y,x*y])) == (2,2,2)
        @test polyerror(DACE.hessian(x*x+y)[1,1], 2) == 0
        @test polyerror(DACE.jacobian([x*x+y,x*y])[1,2], 1) == 0
    end

    @testset "Concurrent callers" begin
        DACE.init(4,2)
        x,y = DACE.identity()
        expected = DACE.getCoefficient(sin(1+x)*exp(y), [2,1])
        jobs = [Threads.@spawn begin
            for _ in 1:100
                value = DACE.getCoefficient(sin(1+x)*exp(y), [2,1])
                value == expected || error("Cross-task engine state corruption")
            end
            true
        end for _ in 1:8]
        @test all(fetch, jobs)
    end
    # Check our methods against Julia's numeric API. Combinations of distinct AD
    # scalar types (e.g. DA with ForwardDiff.Dual) need their own promotion policy.
    standard(m) = m.module in (DACE, Base, Core, SpecialFunctions) || parentmodule(m.module) === Base
    @test isempty(filter(pair -> all(standard, pair), Test.detect_ambiguities(DACE; recursive=true)))
end
