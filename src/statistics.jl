"All multi-indices of total degree <= order, ordered by degree then reverse lexicographically."
function getMultiIndices(no::Integer, nv::Integer)
    no >= 0 && nv >= 1 || throw(ArgumentError("Invalid order or dimension"))
    no <= typemax(UInt32) && nv <= typemax(Int) || throw(ArgumentError("Dimensions exceed index range"))
    result = Vector{UInt32}[]
    current = zeros(UInt32, nv)
    function visit(i, remaining)
        if i == nv
            current[i] = remaining
            push!(result, copy(current))
        else
            for value in remaining:-1:0
                current[i] = value
                visit(i+1, remaining-value)
            end
        end
    end
    for degree in 0:no
        visit(1, degree)
    end
    result
end
function getRawMoments(mgf::DA, no::Integer)
    lock(engine_lock) do
        valid(mgf)
        0 <= no <= getMaxOrder() || throw(ArgumentError("Moment order exceeds initialized order"))
        indices = getMultiIndices(no, getMaxVariables())
        moments = [getCoefficient(mgf, jj)*prod(j -> Float64(factorial(big(j))), jj) for jj in indices]
        indices, moments
    end
end
function getCentralMoments(mgf::DA, no::Integer)
    lock(engine_lock) do
        valid(mgf)
        mean = linear(mgf)
        shift = sum(mean[i]*variable(i) for i in eachindex(mean))
        getRawMoments(exp(-shift)*mgf, no)
    end
end
