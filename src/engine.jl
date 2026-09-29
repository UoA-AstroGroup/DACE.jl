const lib = normpath(joinpath(@__DIR__, "../deps/usr/lib/libdacecore.$(Libdl.dlext)"))
const engine_lock = ReentrantLock()
const epoch = Ref(UInt(0))
const initialized = Ref(false)
const max_order = Ref(0)
const variables = Ref(0)
const monomials = Ref(0)
const truncation_stack = UInt32[]
const monomial_words = Ref(0)

"Layout of DACEDA with DACE_MEMORY_DYNAMIC; verified against the compiled header."
struct C_DA
    len::Cuint
    capacity::Cuint
    mem::Ptr{Cvoid}
end

# Julia owns the storage passed to DACE's C engine.
mutable struct DA <: Real
    storage::Base.RefValue{C_DA}
    buffer::Vector{UInt64}
    generation::UInt
    alive::Bool
end

# Passing the owner to ccall preserves both its descriptor and coefficient buffer.
Base.cconvert(::Type{Ptr{C_DA}}, a::DA) = a
function Base.unsafe_convert(::Type{Ptr{C_DA}}, a::DA)
    d = a.storage[]
    a.storage[] = C_DA(d.len, d.capacity, pointer(a.buffer))
    Base.unsafe_convert(Ptr{C_DA}, a.storage)
end

struct DACEError <: Exception
    code::UInt32
    message::String
end
Base.showerror(io::IO, e::DACEError) = print(io, "DACE error ", e.code, ": ", e.message)

function __init__()
    isfile(lib) || error("DACE's C engine is missing. Run Pkg.build(\"DACE\") first.")
    ccall((:daceDirectMemoryModel, lib), Cint, ()) == 3 || error("Expected DYNAMIC memory model")
    ccall((:daceDirectSizeofDA, lib), Csize_t, ()) == sizeof(C_DA) || error("DACEDA size mismatch")
    ccall((:daceDirectOffsetMem, lib), Csize_t, ()) == fieldoffset(C_DA, 3) || error("DACEDA layout mismatch")
    bytes = ccall((:daceDirectSizeofMonomial, lib), Csize_t, ())
    alignment = ccall((:daceDirectAlignMonomial, lib), Csize_t, ())
    bytes % sizeof(UInt64) == 0 && alignment <= sizeof(UInt64) || error("Unsupported monomial layout")
    monomial_words[] = bytes ÷ sizeof(UInt64)
end

function check_error()
    code = ccall((:daceGetError, lib), Cuint, ())
    code == 0 && return
    message = unsafe_string(ccall((:daceGetErrorMSG, lib), Cstring, ()))
    severity = ccall((:daceGetErrorX, lib), Cuint, ())
    ccall((:daceClearError, lib), Cvoid, ())
    if severity < 6
        @warn message code
        return
    end
    if severity >= 9
        initialized[] = false
        epoch[] += 1
    end
    throw(DACEError(code, message))
end

ready() = initialized[] || throw(ArgumentError("Call init(order, variables) first"))
function valid(a::DA)
    ready()
    a.alive || throw(ArgumentError("Polynomial has been closed"))
    a.generation == epoch[] || throw(ArgumentError("Polynomial belongs to a previous initialization"))
    return a
end

function init(no::Integer, nv::Integer)
    no >= 1 && nv >= 1 || throw(ArgumentError("Order and variable count must be positive"))
    no <= typemax(Cuint) && nv <= typemax(Cuint) || throw(ArgumentError("Order or variable count exceeds the C API"))
    # The C engine indexes these tables with unsigned int. Bound the exponent
    # first so malformed inputs cannot allocate a huge BigInt during validation.
    cld(nv, 2)*log2(Float64(no)+1) < 32 || throw(ArgumentError("DACE lookup table exceeds its index range"))
    binomial(big(no)+nv, nv) <= typemax(Cuint) || throw(ArgumentError("Too many monomials for DACE"))
    lock(engine_lock) do
        initialized[] = false
        epoch[] += 1
        ccall((:daceInitialize, lib), Cvoid, (Cuint, Cuint), no, nv)
        check_error()
        max_order[], variables[] = no, nv
        monomials[] = ccall((:daceGetMaxMonomials, lib), Cuint, ())
        initialized[] = true
        empty!(truncation_stack)
    end
    return nothing
end

function release!(a::DA)
    if a.alive
        a.alive = false
        a.storage[] = C_DA(0, 0, C_NULL)
        a.buffer = UInt64[]
    end
    return nothing
end
Base.close(a::DA) = lock(() -> release!(a), engine_lock)

function allocate(capacity::Integer=monomials[])
    ready()
    # DYNAMIC never reallocates a caller's buffer. Julia accounts for its full
    # size and reclaims it without a native finalizer. Never call daceFreeDA on it.
    capacity = Int(clamp(capacity, 1, monomials[]))
    buffer = Vector{UInt64}(undef, Base.checked_mul(capacity, monomial_words[]))
    DA(Ref(C_DA(0, capacity, pointer(buffer))), buffer, epoch[], true)
end

function reserve!(a::DA, capacity::Integer)
    d = a.storage[]
    capacity <= d.capacity && return a
    capacity = min(monomials[], max(capacity, 2Int(d.capacity)))
    resize!(a.buffer, Base.checked_mul(capacity, monomial_words[]))
    a.storage[] = C_DA(d.len, capacity, pointer(a.buffer))
    a
end

DA() = DA(0.0)
function _constant(value::Real, capacity::Integer=1)
    lock(engine_lock) do
        a = allocate(capacity)
        ccall((:daceCreateConstant, lib), Cvoid, (Ptr{C_DA}, Cdouble), a, value)
        check_error()
        a
    end
end
DA(value::Real) = _constant(value)
DA(a::DA) = copy(a)
function DA(index::Integer, coefficient::Real)
    lock(engine_lock) do
        ready()
        0 <= index <= variables[] || throw(ArgumentError("Variable index out of bounds"))
        a = allocate(1)
        ccall((:daceCreateVariable, lib), Cvoid, (Ptr{C_DA}, Cuint, Cdouble), a, index, coefficient)
        check_error()
        a
    end
end
function variable(i::Integer)
    i >= 1 || throw(ArgumentError("Variable indices start at one"))
    DA(i, 1.0)
end

function cons(a::DA)
    lock(engine_lock) do
        valid(a)
        result = ccall((:daceGetConstant, lib), Cdouble, (Ptr{C_DA},), a)
        check_error()
        result
    end
end
cons(a::Real) = a

function linear(a::DA)
    lock(engine_lock) do
        valid(a)
        result = zeros(variables[])
        ccall((:daceGetLinear, lib), Cvoid, (Ptr{C_DA}, Ptr{Cdouble}), a, result)
        check_error()
        result
    end
end
linear(v::AbstractVector{<:DA}) = reduce(vcat, permutedims.(linear.(v)))

function coefficient(a::DA, exponents::AbstractVector{<:Integer})
    lock(engine_lock) do
        valid(a)
        length(exponents) == variables[] || throw(DimensionMismatch("One exponent per variable is required"))
        all(>=(0), exponents) && sum(big, exponents) <= max_order[] || throw(ArgumentError("Invalid exponents"))
        jj = Cuint.(exponents)
        result = ccall((:daceGetCoefficient, lib), Cdouble, (Ptr{C_DA}, Ptr{Cuint}), a, jj)
        check_error()
        result
    end
end
