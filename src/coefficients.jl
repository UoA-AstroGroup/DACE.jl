getMaxMonomials() = monomials[]
for (name, c_name, T) in ((:getEps, :daceGetEpsilon, Cdouble),
        (:getEpsMac, :daceGetMachineEpsilon, Cdouble), (:getTO, :daceGetTruncationOrder, Cuint))
    @eval function $name()
        lock(engine_lock) do
            ready()
            ccall(($(QuoteNode(c_name)), lib), $T, ())
        end
    end
end
function setTO(n::Integer)
    lock(engine_lock) do
        ready()
        0 <= n <= max_order[] || throw(ArgumentError("Truncation order must lie between zero and the initialized order"))
        old = ccall((:daceSetTruncationOrder, lib), Cuint, (Cuint,), n)
        check_error()
        Int(old)
    end
end
function pushTO(n::Integer)
    lock(engine_lock) do
        old = setTO(n)
        push!(truncation_stack, old)
    end
    nothing
end
function popTO()
    lock(engine_lock) do
        isempty(truncation_stack) && throw(ArgumentError("Truncation-order stack is empty"))
        setTO(pop!(truncation_stack))
    end
    nothing
end

"A coefficient and its nonnegative exponent vector, independent of engine lifetime."
struct Monomial
    coefficient::Float64
    exponents::Vector{UInt32}
end
getCoefficient(m::Monomial) = m.coefficient
getExponents(m::Monomial) = copy(m.exponents)
order(m::Monomial) = sum(m.exponents)
function exponent_vector(jj::AbstractVector{<:Integer})
    all(j -> 0 <= j <= typemax(Cuint), jj) || throw(ArgumentError("Invalid monomial exponents"))
    out = zeros(Cuint, variables[])
    for (i,j) in enumerate(jj)
        i > length(out) && break
        out[i] = j
    end
    out
end
function getCoefficient(a::DA, jj::AbstractVector{<:Integer})
    lock(engine_lock) do
        valid(a)
        exponents = exponent_vector(jj)
        sum(UInt64, exponents) > max_order[] && return 0.0
        value = ccall((:daceGetCoefficient, lib), Cdouble, (Ptr{C_DA}, Ptr{Cuint}), a, exponents)
        check_error()
        value
    end
end
function setCoefficient!(a::DA, jj::AbstractVector{<:Integer}, value::Real)
    lock(engine_lock) do
        valid(a)
        exponents = exponent_vector(jj)
        sum(UInt64, exponents) <= max_order[] || throw(ArgumentError("Monomial exceeds initialized order"))
        reserve!(a, Int(a.storage[].len)+1)
        ccall((:daceSetCoefficient, lib), Cvoid, (Ptr{C_DA}, Ptr{Cuint}, Cdouble), a, exponents, value)
        check_error()
        a
    end
end
function getMonomial(a::DA, pos::Integer)
    lock(engine_lock) do
        valid(a)
        1 <= pos <= a.storage[].len || throw(BoundsError(a, pos))
        jj, c = zeros(Cuint, variables[]), Ref{Cdouble}(0)
        ccall((:daceGetCoefficientAt, lib), Cvoid, (Ptr{C_DA}, Cuint, Ptr{Cuint}, Ref{Cdouble}), a, pos, jj, c)
        check_error()
        Monomial(c[], jj)
    end
end
function getMonomials(a::DA)
    lock(engine_lock) do
        valid(a)
        sort!([getMonomial(a, i) for i in 1:a.storage[].len]; by=order, alg=Base.Sort.MergeSort)
    end
end
function trim(a::DA, low::Integer, high::Integer=getMaxOrder())
    0 <= low <= typemax(Cuint) && 0 <= high <= typemax(Cuint) || throw(ArgumentError("Invalid order range"))
    lock(engine_lock) do
        valid(a)
        low > high && return DA(0.0)
        result = allocate()
        ccall((:daceTrim, lib), Cvoid, (Ptr{C_DA}, Cuint, Cuint, Ptr{C_DA}), a, low, high, result)
        check_error()
        result
    end
end
function multiplyMonomials(a::DA, b::DA)
    lock(engine_lock) do
        valid(a); valid(b)
        # Upstream daceMultiplyMonomials does not set the output length and scans
        # capacity rather than used terms. Use the public coefficient API instead.
        result = DA(0.0)
        for i in 1:a.storage[].len
            monomial = getMonomial(a, i)
            value = monomial.coefficient * getCoefficient(b, monomial.exponents)
            iszero(value) || setCoefficient!(result, monomial.exponents, value)
        end
        result
    end
end
function divide(a::DA, i::Integer, p::Integer=1)
    lock(engine_lock) do
        valid(a)
        1 <= i <= variables[] && 0 <= p <= max_order[] || throw(ArgumentError("Invalid variable power"))
        result = allocate()
        ccall((:daceDivideByVariable, lib), Cvoid, (Ptr{C_DA}, Cuint, Cuint, Ptr{C_DA}), a, i, p, result)
        check_error()
        result
    end
end
const integ = integrate
for (name, c_name) in ((:deriv, :daceDifferentiate), (:integrate, :daceIntegrate))
    @eval function $name(a::DA, counts::AbstractVector{<:Integer})
        all(>=(0), counts) || throw(ArgumentError("Derivative/integral counts must be nonnegative"))
        lock(engine_lock) do
            result = copy(a)
            for (i,n) in enumerate(counts)
                i > variables[] && break
                for _ in 1:n
                    ccall(($(QuoteNode(c_name)), lib), Cvoid, (Cuint, Ptr{C_DA}, Ptr{C_DA}), i, result, result)
                    check_error()
                end
            end
            result
        end
    end
end
function orderNorm(a::DA, v::Integer=0, p::Integer=0)
    lock(engine_lock) do
        valid(a)
        0 <= v <= variables[] && 0 <= p <= typemax(Cuint) || throw(ArgumentError("Invalid norm arguments"))
        out = zeros(max_order[]+1)
        ccall((:daceOrderedNorm, lib), Cvoid, (Ptr{C_DA}, Cuint, Cuint, Ptr{Cdouble}), a, v, p, out)
        check_error()
        out
    end
end
function estimNorm(a::DA, v::Integer=0, p::Integer=0, n::Integer=getMaxOrder())
    lock(engine_lock) do
        valid(a)
        0 <= v <= variables[] && 0 <= p <= typemax(Cuint) && 0 <= n < typemax(Cuint) || throw(ArgumentError("Invalid norm arguments"))
        out = zeros(n+1)
        ccall((:daceEstimate, lib), Cvoid, (Ptr{C_DA}, Cuint, Cuint, Ptr{Cdouble}, Ptr{Cdouble}, Cuint), a, v, p, out, C_NULL, n)
        check_error()
        out
    end
end
struct Interval
    m_lb::Float64
    m_ub::Float64
end
function bound(a::DA)
    lock(engine_lock) do
        valid(a)
        lo, hi = Ref{Cdouble}(0), Ref{Cdouble}(0)
        ccall((:daceGetBounds, lib), Cvoid, (Ptr{C_DA}, Ref{Cdouble}, Ref{Cdouble}), a, lo, hi)
        check_error()
        Interval(lo[], hi[])
    end
end
function random(filling::Real=-1.0)
    lock(engine_lock) do
        result = allocate()
        ccall((:daceCreateRandom, lib), Cvoid, (Ptr{C_DA}, Cdouble), result, filling)
        check_error()
        result
    end
end
function toString(a::DA)
    lock(engine_lock) do
        valid(a)
        width = Int(ccall((:daceDirectStringLength, lib), Csize_t, ()))
        bytes = zeros(UInt8, width*(Int(a.storage[].len)+2))
        n = Ref{Cuint}(0)
        ccall((:daceWrite, lib), Cvoid, (Ptr{C_DA}, Ptr{UInt8}, Ref{Cuint}), a, bytes, n)
        check_error()
        lines = [String(bytes[(i-1)*width+1:(i-1)*width+findfirst(iszero, @view(bytes[(i-1)*width+1:i*width]))-1]) for i in 1:Int(n[])]
        join(lines, '\n')*"\n"
    end
end
toString(m::Monomial) = string(m.coefficient, " ", m.exponents, '\n')
function Base.show(io::IO, a::DA)
    if !a.alive || a.generation != epoch[]
        print(io, "DA(inactive)")
    elseif get(io, :compact, false)
        print(io, "DA(", cons(a), "; ", a.storage[].len, " terms)")
    else
        print(io, toString(a))
    end
end
Base.show(io::IO, m::Monomial) = print(io, toString(m))
