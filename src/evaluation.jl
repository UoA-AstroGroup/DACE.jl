"C-generated evaluation tree, interpreted by Julia for numbers or polynomials."
struct CompiledMap
    coefficients::Matrix{Float64}
    levels::Vector{Int}
    indices::Vector{Int}
    depth::Int
    nvars::Int
    generation::UInt
end
const compiledDA = CompiledMap
CompiledMap(a::Union{DA,AbstractVector{<:DA}}) = compile(a)
getDim(map::CompiledMap) = size(map.coefficients, 1)
getOrd(map::CompiledMap) = map.depth
getVars(map::CompiledMap) = map.nvars
getTerms(map::CompiledMap) = length(map.levels)

function compile(polynomials::AbstractVector{<:DA})
    isempty(polynomials) && throw(ArgumentError("Cannot compile an empty map"))
    lock(engine_lock) do
        # Lazy arrays may create a fresh DA at each getindex. Retain the actual
        # values, not just the input array, while C traverses the pointer table.
        owners = polynomials isa Vector{DA} ? polynomials : collect(polynomials)
        foreach(valid, owners)
        n = length(owners)
        data = Vector{Float64}(undef, (n + 2) * monomials[])
        terms, nv, no = Ref{Cuint}(0), Ref{Cuint}(0), Ref{Cuint}(0)
        # The pointer array alone does not own the polynomials' Julia buffers.
        GC.@preserve owners begin
            pointers = [Base.unsafe_convert(Ptr{C_DA}, p) for p in owners]
            ccall((:daceEvalTree, lib), Cvoid,
                (Ptr{Ptr{C_DA}}, Cuint, Ptr{Cdouble}, Ref{Cuint}, Ref{Cuint}, Ref{Cuint}),
                pointers, n, data, terms, nv, no)
        end
        check_error()
        resize!(data, (n + 2) * Int(terms[]))
        tree = reshape(data, n + 2, :)
        # DACE encodes traversal indices as doubles. Decode once, not at every
        # evaluation: checked Float64-to-Int conversions dominate small maps.
        CompiledMap(tree[3:end, :], Int.(tree[1, :]), Int.(tree[2, :]),
            Int(no[]), Int(nv[]), epoch[])
    end
end
compile(a::DA) = compile([a])

function evaluate!(result::Vector{Float64}, map::CompiledMap, args::AbstractVector{<:Real}, work::Vector{Float64})
    Base.require_one_based_indexing(args)
    any(x -> x isa DA, args) && throw(ArgumentError("Use evaluate for polynomial composition"))
    (Base.mightalias(result, args) || Base.mightalias(work, args) || Base.mightalias(result, work)) &&
        throw(ArgumentError("Output, inputs and workspace must not alias"))
    # DACE treats unspecified variables as zero and ignores extra coordinates.
    if length(args) < map.nvars
        return evaluate!(result, map, vcat(args, zeros(map.nvars-length(args))), work)
    end
    n = size(map.coefficients, 1)
    length(result) == n && length(work) >= map.depth + 1 || throw(DimensionMismatch("Invalid workspace"))
    work[1] = 1.0
    @inbounds for j in 1:n
        result[j] = map.coefficients[j, 1]
    end
    @inbounds for i in 2:length(map.levels)
        depth, var = map.levels[i], map.indices[i]
        work[depth + 1] = work[depth] * args[var]
        for j in 1:n
            result[j] += work[depth + 1] * map.coefficients[j, i]
        end
    end
    result
end

function evaluate(map::CompiledMap, args::AbstractVector{<:Real})
    any(x -> x isa DA, args) && return evaluate(map, DA.(args))
    evaluate!(zeros(size(map.coefficients, 1)), map, args, zeros(map.depth + 1))
end

function evaluate(map::CompiledMap, args::AbstractVector{<:DA})
    Base.require_one_based_indexing(args)
    lock(engine_lock) do
        foreach(valid, args)
        map.generation == epoch[] || throw(ArgumentError("Compiled map belongs to a previous initialization"))
        if length(args) < map.nvars
            return evaluate(map, vcat(args, [DA(0.0) for _ in 1:map.nvars-length(args)]))
        end
        n = size(map.coefficients, 1)
        result = [_constant(map.coefficients[j, 1], monomials[]) for j in 1:n]
        work = [allocate() for _ in 0:map.depth]
        scratch = allocate()
        try
            ccall((:daceCreateConstant, lib), Cvoid, (Ptr{C_DA}, Cdouble), work[1], 1.0)
            for i in 2:length(map.levels)
                depth, var = map.levels[i], map.indices[i]
                ccall((:daceMultiply, lib), Cvoid, (Ptr{C_DA}, Ptr{C_DA}, Ptr{C_DA}),
                    work[depth], args[var], work[depth + 1])
                check_error()
                for j in 1:n
                    c = map.coefficients[j, i]
                    iszero(c) && continue
                    # daceWeightedSum forbids aliasing either input with its output.
                    ccall((:daceWeightedSum, lib), Cvoid,
                        (Ptr{C_DA}, Cdouble, Ptr{C_DA}, Cdouble, Ptr{C_DA}),
                        result[j], 1.0, work[depth + 1], c, scratch)
                    check_error()
                    result[j], scratch = scratch, result[j]
                end
            end
        finally
            foreach(close, work)
            close(scratch)
        end
        result
    end
end
evaluate(a::DA, args::AbstractVector{<:Real}) = only(evaluate(compile(a), args))
evaluate(a::AbstractVector{<:DA}, args::AbstractVector{<:Real}) = evaluate(compile(a), args)
evalScalar(a::DA, value::Real) = evaluate(a, [value])
evalScalar(a::Union{CompiledMap,AbstractVector{<:DA}}, value::Real) = evaluate(a, [value])
function evaluate(map::CompiledMap, args::AbstractVector{<:Real}, result::AbstractVector)
    values = evaluate(map, args)
    length(values) == length(result) || throw(DimensionMismatch("Output has incorrect length"))
    copyto!(result, values)
    result
end
evaluate(map::CompiledMap, args::AbstractVector{<:AbstractFloat}, result::Vector{Float64}) =
    evaluate!(result, map, args, zeros(map.depth + 1))

function _linear_transform(A::AbstractMatrix{<:Real}, x::AbstractVector{<:DA})
    size(A,2) == length(x) || throw(DimensionMismatch("Linear map dimensions differ"))
    result = [allocate() for _ in axes(A,1)]
    scratch = allocate()
    for j in eachindex(x), i in eachindex(result)
        c = A[i,j]
        iszero(c) && continue
        ccall((:daceWeightedSum, lib), Cvoid,
            (Ptr{C_DA}, Cdouble, Ptr{C_DA}, Cdouble, Ptr{C_DA}), result[i], 1.0, x[j], c, scratch)
        check_error()
        result[i], scratch = scratch, result[i]
    end
    result
end

"Invert a square polynomial map through the current truncation order."
function invert(f::AbstractVector{<:DA})
    lock(engine_lock) do
        ready()
        foreach(valid, f)
        1 <= length(f) <= variables[] || throw(DimensionMismatch("Map dimension exceeds the independent variables"))
        if length(f) < variables[]
            # Complete a partial map with identity coordinates for its parameters.
            return invert(vcat(f, variable.(length(f)+1:variables[])))[1:length(f)]
        end
        saved = ccall((:daceGetTruncationOrder, lib), Cuint, ())
        x = variable.(1:variables[])
        c = cons.(f)
        A = linear(f)
        Ainv = lu(A) \ Matrix{Float64}(I, length(f), length(f))
        linear_inverse = _linear_transform(Ainv, x)
        nonlinear = compile(_linear_transform(Ainv, trim.(f, 2)))
        result = linear_inverse
        try
            # For f(x)=c+A*x+h(x), solve g=A^-1*y-A^-1*h(g) one degree at a time.
            for degree in 2:Int(saved)
                ccall((:daceSetTruncationOrder, lib), Cuint, (Cuint,), degree)
                result = linear_inverse - evaluate(nonlinear, result)
            end
            result = evaluate(compile(result), x - c)
        finally
            ccall((:daceSetTruncationOrder, lib), Cuint, (Cuint,), saved)
        end
        check_error()
        result
    end
end
