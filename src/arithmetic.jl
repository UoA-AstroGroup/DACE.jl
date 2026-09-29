for (julia_name, c_name) in ((:copy, :daceCopy), (:sin, :daceSine), (:cos, :daceCosine),
        (:tan, :daceTangent), (:asin, :daceArcSine), (:acos, :daceArcCosine), (:atan, :daceArcTangent),
        (:exp, :daceExponential), (:log, :daceLogarithm),
        (:sinh, :daceHyperbolicSine), (:cosh, :daceHyperbolicCosine), (:tanh, :daceHyperbolicTangent),
        (:asinh, :daceHyperbolicArcSine), (:acosh, :daceHyperbolicArcCosine), (:atanh, :daceHyperbolicArcTangent),
        (:log2, :daceLogarithm2), (:log10, :daceLogarithm10),
        (:round, :daceRound), (:trunc, :daceTruncate), (:inv, :daceMultiplicativeInverse))
    @eval function Base.$julia_name(a::DA)
        lock(engine_lock) do
            valid(a)
            result = allocate($(julia_name in (:copy, :round, :trunc) ? :(a.storage[].len) : :(monomials[])))
            ccall(($(QuoteNode(c_name)), lib), Cvoid, (Ptr{C_DA}, Ptr{C_DA}), a, result)
            check_error()
            result
        end
    end
end

function Base.sqrt(a::DA)
    lock(engine_lock) do
        valid(a)
        # Upstream rejects every zero constant term; the identically zero
        # polynomial has an exact square root and needs no Taylor expansion.
        a.storage[].len == 0 && return DA(0.0)
        result = allocate()
        ccall((:daceSquareRoot, lib), Cvoid, (Ptr{C_DA}, Ptr{C_DA}), a, result)
        check_error()
        result
    end
end

for (op, c_name) in ((:+, :daceAdd), (:-, :daceSubtract), (:*, :daceMultiply), (:/, :daceDivide))
    @eval function Base.$op(a::DA, b::DA)
        lock(engine_lock) do
            valid(a); valid(b)
            capacity = $(op in (:+, :-) ? :(Int(a.storage[].len)+Int(b.storage[].len)) :
                         op == :* ? :(min(UInt64(a.storage[].len)*b.storage[].len, UInt64(monomials[]))) :
                         :(b.storage[].len == 1 && !iszero(cons(b)) ? a.storage[].len : monomials[]))
            result = allocate(capacity)
            ccall(($(QuoteNode(c_name)), lib), Cvoid, (Ptr{C_DA}, Ptr{C_DA}, Ptr{C_DA}), a, b, result)
            check_error()
            result
        end
    end
end

for (op, c_name, reverse_name) in ((:+, :daceAddDouble, :daceAddDouble),
        (:-, :daceSubtractDouble, :daceDoubleSubtract), (:*, :daceMultiplyDouble, :daceMultiplyDouble),
        (:/, :daceDivideDouble, :daceDoubleDivide))
    @eval begin
        function Base.$op(a::DA, b::Real)
            lock(engine_lock) do
                valid(a)
                result = allocate($(op in (:+, :-) ? :(Int(a.storage[].len)+1) : :(a.storage[].len)))
                ccall(($(QuoteNode(c_name)), lib), Cvoid, (Ptr{C_DA}, Cdouble, Ptr{C_DA}), a, b, result)
                check_error()
                result
            end
        end
        function Base.$op(b::Real, a::DA)
            lock(engine_lock) do
                valid(a)
                result = allocate($(op in (:+, :-) ? :(Int(a.storage[].len)+1) : op == :/ ? :(monomials[]) : :(a.storage[].len)))
                ccall(($(QuoteNode(reverse_name)), lib), Cvoid, (Ptr{C_DA}, Cdouble, Ptr{C_DA}), a, b, result)
                check_error()
                result
            end
        end
    end
end

function Base.:^(a::DA, n::Integer)
    typemin(Cint) < n <= typemax(Cint) || throw(ArgumentError("Exponent exceeds supported C int range"))
    lock(engine_lock) do
        valid(a)
        result = allocate()
        ccall((:dacePower, lib), Cvoid, (Ptr{C_DA}, Cint, Ptr{C_DA}), a, n, result)
        check_error()
        result
    end
end
function Base.:^(a::DA, n::AbstractFloat)
    lock(engine_lock) do
        valid(a)
        result = allocate()
        ccall((:dacePowerDouble, lib), Cvoid, (Ptr{C_DA}, Cdouble, Ptr{C_DA}), a, n, result)
        check_error()
        result
    end
end
Base.:-(a::DA) = -1.0 * a
Base.:+(a::DA) = a
Base.zero(::Type{DA}) = DA(0.0)
Base.one(::Type{DA}) = DA(1.0)
Base.zero(::DA) = DA(0.0)
Base.one(::DA) = DA(1.0)
Base.convert(::Type{DA}, x::Real) = DA(x)
Base.convert(::Type{DA}, x::DA) = x
Base.promote_rule(::Type{DA}, ::Type{T}) where {T<:Real} = DA
Base.conj(a::DA) = a
Base.real(a::DA) = a
Base.abs2(a::DA) = a * a
function Base.iszero(a::DA)
    lock(engine_lock) do
        valid(a)
        # The C engine removes zero coefficients. A coefficient norm cannot
        # identify zero reliably when coefficients contain NaN.
        a.storage[].len == 0
    end
end
Base.:(==)(a::DA, b::DA) = cons(a) == cons(b)
Base.:(==)(a::DA, b::Real) = cons(a) == b
Base.:(==)(a::Real, b::DA) = b == a
Base.:(==)(a::DA, b::AbstractIrrational) = cons(a) == b
Base.:(==)(a::AbstractIrrational, b::DA) = a == cons(b)
Base.abs(a::DA) = signbit(cons(a)) ? -a : copy(a)
Base.signbit(a::DA) = signbit(cons(a))
Base.sign(a::DA) = sign(cons(a))
Base.float(a::DA) = a
Base.eps(a::DA) = eps(cons(a))
Base.eps(::Type{DA}) = eps(Float64)
Base.hash(a::DA, h::UInt) = hash(cons(a), h)
Base.Float64(a::DA) = cons(a)
for op in (:<, :<=, :>, :>=, :isless, :isequal)
    @eval begin
        Base.$op(a::DA, b::DA) = Base.$op(cons(a), cons(b))
        Base.$op(a::DA, b::Real) = Base.$op(cons(a), b)
        Base.$op(a::Real, b::DA) = Base.$op(a, cons(b))
    end
end

# Resolve the intersections with Base's special floating-point ordering methods.
for op in (:isless, :isequal)
    @eval begin
        Base.$op(a::DA, b::AbstractFloat) = Base.$op(cons(a), b)
        Base.$op(a::AbstractFloat, b::DA) = Base.$op(a, cons(b))
    end
end

for (name, c_name) in ((:deriv, :daceDifferentiate), (:integrate, :daceIntegrate))
    @eval function $name(a::DA, index::Integer)
        lock(engine_lock) do
            valid(a)
            1 <= index <= variables[] || throw(ArgumentError("Variable index out of bounds"))
            result = allocate(a.storage[].len)
            ccall(($(QuoteNode(c_name)), lib), Cvoid, (Cuint, Ptr{C_DA}, Ptr{C_DA}), index, a, result)
            check_error()
            result
        end
    end
end

function norm(a::DA, p::Real=0)
    isinteger(p) && 0 <= p <= typemax(Cuint) || throw(ArgumentError("Use a nonnegative integer coefficient norm"))
    lock(engine_lock) do
        valid(a)
        result = ccall((:daceNorm, lib), Cdouble, (Ptr{C_DA}, Cuint), a, p)
        check_error()
        result
    end
end
function setEps(value::Real)
    isfinite(value) && value >= 0 || throw(ArgumentError("Epsilon must be finite and nonnegative"))
    lock(engine_lock) do
        ready()
        previous = ccall((:daceSetEpsilon, lib), Cdouble, (Cdouble,), value)
        check_error()
        previous
    end
end
getMaxOrder() = max_order[]
getMaxVariables() = variables[]

powi(a::DA, p::Integer) = a^p
powd(a::DA, p::Real) = a^Float64(p)
Base.:^(a::DA, p::Rational) = a^Float64(p)
# Keep the full exponent polynomial, rather than dropping its nonconstant part.
Base.:^(a::Real, b::DA) = exp(log(a)*b)
Base.:^(a::DA, b::DA) = exp(log(a)*b)
Base.:^(::Irrational{:ℯ}, b::DA) = exp(b)
Base.log(::Irrational{:ℯ}, a::DA) = log(a)

for (name, c_name) in ((:isnan, :daceIsNan), (:isinf, :daceIsInf))
    @eval function Base.$name(a::DA)
        lock(engine_lock) do
            valid(a)
            result = ccall(($(QuoteNode(c_name)), lib), Cuint, (Ptr{C_DA},), a)
            check_error()
            result != 0
        end
    end
end
Base.isfinite(a::DA) = !isnan(a) && !isinf(a)

for (name, c_name) in ((:atan, :daceArcTangent2), (:hypot, :daceHypotenuse))
    @eval function Base.$name(a::DA, b::DA)
        lock(engine_lock) do
            valid(a); valid(b)
            $(name == :hypot) && a.storage[].len == 0 && b.storage[].len == 0 && return DA(0.0)
            result = allocate()
            ccall(($(QuoteNode(c_name)), lib), Cvoid, (Ptr{C_DA}, Ptr{C_DA}, Ptr{C_DA}), a, b, result)
            check_error()
            result
        end
    end
    @eval Base.$name(a::DA, b::Real) = Base.$name(a, DA(b))
    @eval Base.$name(a::Real, b::DA) = Base.$name(DA(a), b)
end

function root(a::DA, p::Integer=2)
    typemin(Cint) < p <= typemax(Cint) || throw(ArgumentError("Root exceeds C int"))
    p != 0 || throw(DomainError(p, "Zeroth root is undefined"))
    lock(engine_lock) do
        valid(a)
        p > 0 && iszero(a) && return DA(0.0)
        result = allocate()
        ccall((:daceRoot, lib), Cvoid, (Ptr{C_DA}, Cint, Ptr{C_DA}), a, p, result)
        check_error()
        result
    end
end
isrt(a::DA) = root(a, -2)
icrt(a::DA) = root(a, -3)
Base.cbrt(a::DA) = root(a, 3)
sqr(a::DA) = a*a
function Base.log(b::Real, a::DA)
    lock(engine_lock) do
        valid(a)
        result = allocate()
        ccall((:daceLogarithmBase, lib), Cvoid, (Ptr{C_DA}, Cdouble, Ptr{C_DA}), a, b, result)
        check_error()
        result
    end
end
function Base.mod(a::DA, p::Real)
    lock(engine_lock) do
        valid(a)
        result = allocate()
        ccall((:daceModulo, lib), Cvoid, (Ptr{C_DA}, Cdouble, Ptr{C_DA}), a, p, result)
        check_error()
        result
    end
end

for (name, c_name) in ((:erf, :daceErrorFunction), (:erfc, :daceComplementaryErrorFunction),
        (:gamma, :daceGammaFunction), (:loggamma, :daceLogGammaFunction))
    @eval const $name = SpecialFunctions.$name
    @eval function SpecialFunctions.$name(a::DA)
        lock(engine_lock) do
            valid(a)
            result = allocate()
            ccall(($(QuoteNode(c_name)), lib), Cvoid, (Ptr{C_DA}, Ptr{C_DA}), a, result)
            check_error()
            result
        end
    end
end
for (name, c_name) in ((:besselj, :daceBesselJFunction), (:bessely, :daceBesselYFunction))
    @eval const $name = SpecialFunctions.$name
    @eval function SpecialFunctions.$name(n::Integer, a::DA)
        typemin(Cint) < n <= typemax(Cint) || throw(ArgumentError("Bessel order exceeds C int"))
        lock(engine_lock) do
            valid(a)
            result = allocate()
            ccall(($(QuoteNode(c_name)), lib), Cvoid, (Ptr{C_DA}, Cint, Ptr{C_DA}), a, n, result)
            check_error()
            result
        end
    end
end
for (name, scaledname, c_name) in ((:besseli, :besselix, :daceBesselIFunction), (:besselk, :besselkx, :daceBesselKFunction))
    @eval const $name = SpecialFunctions.$name
    for (fn, scaled) in ((name, false), (scaledname, true))
        @eval function SpecialFunctions.$fn(n::Integer, a::DA)
            typemin(Cint) < n <= typemax(Cint) || throw(ArgumentError("Bessel order exceeds C int"))
            lock(engine_lock) do
                valid(a)
                result = allocate()
                ccall(($(QuoteNode(c_name)), lib), Cvoid, (Ptr{C_DA}, Cint, Bool, Ptr{C_DA}), a, n, $scaled, result)
                check_error()
                result
            end
        end
    end
end
function PsiFunction(a::DA, n::Integer)
    0 <= n <= typemax(Cuint) || throw(ArgumentError("Invalid polygamma order"))
    lock(engine_lock) do
        valid(a)
        result = allocate()
        ccall((:dacePsiFunction, lib), Cvoid, (Ptr{C_DA}, Cuint, Ptr{C_DA}), a, n, result)
        check_error()
        result
    end
end
SpecialFunctions.polygamma(n::Integer, a::DA) = PsiFunction(a, n)
SpecialFunctions.digamma(a::DA) = PsiFunction(a, 0)
