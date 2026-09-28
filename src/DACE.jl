baremodule DACE
    using Base

    # Julia 1.12+ gives ordinary modules an `eval` binding of type Core.EvalInto.
    # CxxWrap needs a generic function to add DACE's polynomial-evaluation
    # methods. Define the usual module helpers explicitly to preserve that API.
    eval(x) = Core.eval(DACE, x)
    include(path) = Base.include(DACE, path)
    include(mapexpr::Function, path) = Base.include(mapexpr, DACE, path)

    using DACE_jll
    using CxxWrap
    using SpecialFunctions
    using DiffEqBase
    using DocStringExtensions

    DocStringExtensions.@template DEFAULT =
        """
        $(TYPEDSIGNATURES)
        $(DOCSTRING)
        """

    mutable struct Interval
        m_lb::Float64
        m_ub::Float64
    end

    # load the C++ interface and initialize it
    @wrapmodule(() -> libdace, :define_julia_module)
    function __init__()
        @initcxx
    end

    # include the file that contains documentation for some of the DACE C++ wrapped functions
    include("docs.jl")

    # ========================= #
    # extend DACE functionality #
    # ========================= #

    # type alias for convenience
    const DAOrF64 = Union{DA, Float64}

    # ------------------- #
    # custom constructors #
    # ------------------- #

    DA(b::Bool) = DA(b ? 1.0 : 0.0)
    DA(x::Rational) = DA(convert(Float64,x))

    # concrete DA types
    DAAllocated() = DA(0.0)
    DAAllocated(x::Real) = DA(x)

    # AlgebraicVector and AlgebraicMatrix constructors from Julia vectors and matrices
    AlgebraicVector(v::AbstractVector{<:DA}) = AlgebraicVector{DA}(v)
    AlgebraicMatrix(m::AbstractMatrix{<:DA}) = AlgebraicMatrix{DA}(m)
    AlgebraicVector(v::AbstractVector{Float64}) = AlgebraicVector{Float64}(v)
    AlgebraicMatrix(m::AbstractMatrix{Float64}) = AlgebraicMatrix{Float64}(m)

    AlgebraicVector{T}(v::AbstractVector{U}) where {T<:DAOrF64, U<:DAOrF64} = begin
        res = AlgebraicVector{T}(length(v))
        res .= v
        return res
    end

    AlgebraicMatrix{T}(m::AbstractMatrix{U}) where {T<:DAOrF64, U<:DAOrF64} = begin
        res = AlgebraicMatrix{T}(size(m)...)
        res .= m
        return res
    end

    # -------------------------- #
    # overload functions in Base #
    # -------------------------- #

    # additive and multiplicative identities
    Base.zero(::Type{DA}) = DA(0.0)
    Base.one(::Type{DA}) = DA(1.0)

    # similar function for AlgebraicVector and AlgebraicMatrix
    Base.similar(v::AlgebraicVector{<:DA}) = AlgebraicVector{DA}(v)
    Base.similar(m::AlgebraicMatrix{<:DA}) = AlgebraicMatrix{DA}(m)
    Base.similar(v::AlgebraicVector{Float64}) = AlgebraicVector{Float64}(v)
    Base.similar(m::AlgebraicMatrix{Float64}) = AlgebraicMatrix{Float64}(m)

    # promotions and conversions to DA
    @inline Base.promote_rule(::Type{T}, ::Type{R}) where {T<:DA, R<:Real} = T
    @inline Base.promote_rule(::Type{R}, ::Type{T}) where {T<:DA, R<:Real} = T
    # FIXME: how to handle these type promotions properly?
    @inline Base.promote_rule(::Type{T}, ::Type{DAAllocated}) where T<:DA = DAAllocated
    @inline Base.promote_rule(::Type{DAAllocated}, ::Type{T}) where T<:DA = DAAllocated

    for R in (AbstractFloat, AbstractIrrational, Integer, Rational)
        @eval Base.convert(::Type{<:DA}, x::$R) = DA(convert(Float64, x))
        @eval DACE.cons(a::$R) = a
    end

    Base.float(a::DA) = a
    Base.eps(a::DA) = eps(cons(a))
    Base.eps(::Type{T}) where {T<:DA} = eps(Float64)

    # power operators
    Base.:^(da::DA, p::Integer) = DACE.powi(da, p)
    Base.:^(da::DA, p::Float64) = DACE.powd(da, p)

    for R in (AbstractFloat, AbstractIrrational, Integer, Rational)
        @eval Base.:^(a::$R, b::DA) = a^DACE.cons(b)
    end

    Base.isfinite(a::DA) = !(isinf(a) || isnan(a))

    # assignment of vector and matrix elements
    Base.setindex!(v::AlgebraicVector{<:DA}, x::Real, i::Integer) = v[i] = convert(DA, x)
    Base.setindex!(m::AlgebraicMatrix{<:DA}, x::Real, i::Integer, j::Integer) = m[i,j] = convert(DA, x)

    # custom show functions for DA-related objects
    Base.show(io::IO, m::Monomial) = print(io, toString(m))
    Base.show(io::IO, da::DA) = print(io, toString(da))
    Base.show(io::IO, vec::AlgebraicVector{T}) where T<:DAOrF64 = print(io, toString(vec))
    Base.show(io::IO, mat::AlgebraicMatrix{T}) where T<:DAOrF64 = print(io, toString(mat))

    # -------------------------------------- #
    # overload functions in SpecialFunctions #
    # -------------------------------------- #

    # error and complementary error functions
    SpecialFunctions.erf(a::DA) = DACE.erf(a)
    SpecialFunctions.erfc(a::DA) = DACE.erfc(a)

    # -------------------------------- #
    # overload functions in DiffEqBase #
    # -------------------------------- #

    DiffEqBase.value(a::DA) = DACE.cons(a)

    # ---------------------------------------------- #
    # wrappers to handle different subtypes of Array #
    # ---------------------------------------------- #

    # map inversion
    invert(v::AbstractVector{<:DA}) = invert(AlgebraicVector(v))

    # linear part, Jacobian and Hessian
    linear(v::AbstractVector{<:DA}) = linear(AlgebraicVector(v))
    jacobian(v::AbstractVector{<:DA}) = jacobian(AlgebraicVector(v))
    hessian(v::AbstractVector{<:DA}) = hessian(AlgebraicVector(v)) # vector of Hessian matrices
    hess_stack(a::AbstractVector{<:DA}) = stack(hessian(a), dims=3)

    # compilation and evaluation of DA objects
    compile(v::AbstractVector{<:DA}) = compile(AlgebraicVector(v))

    eval(cda::compiledDA, v::AbstractVector{<:DA}) = eval(cda, AlgebraicVector(v))
    eval(a::AbstractVector{<:DA}, v::AbstractVector{<:DA}) = eval(AlgebraicVector(a), AlgebraicVector(v))
    eval(a::AlgebraicVector{<:DA}, v::AbstractVector{<:DA}) = eval(a, AlgebraicVector(v))

    eval(cda::compiledDA, v::AbstractVector{Float64}) = eval(cda, AlgebraicVector(v))
    eval(a::AbstractVector{<:DA}, v::AbstractVector{Float64}) = eval(AlgebraicVector(a), AlgebraicVector(v))
    eval(a::AlgebraicVector{<:DA}, v::AbstractVector{Float64}) = eval(a, AlgebraicVector(v))

    # ------- #
    # exports #
    # ------- #

    export DA, AlgebraicVector, AlgebraicMatrix, compiledDA, Monomial

end  # module DACE
