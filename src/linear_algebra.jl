gradient(a::DA) = [deriv(a, i) for i in 1:getMaxVariables()]
jacobian(a::DA) = permutedims(gradient(a))
jacobian(v::AbstractVector{<:DA}) = [deriv(v[i], j) for i in eachindex(v), j in 1:getMaxVariables()]
hessian(a::DA) = [deriv(deriv(a, i), j) for i in 1:getMaxVariables(), j in 1:getMaxVariables()]
hessian(v::AbstractVector{<:DA}) = hessian.(v)
hess_stack(v::AbstractVector{<:DA}) = stack(hessian(v); dims=3)
LinearAlgebra.det(a::AlgebraicMatrix) = LinearAlgebra.det(a.data)
const det = LinearAlgebra.det
Base.inv(a::AlgebraicMatrix) = AlgebraicMatrix(inv(a.data))
frobenius(a::AbstractMatrix{<:Real}) = sqrt(sum(abs2, a))

"""
    eigh(A)

Eigenvalues and orthonormal eigenvectors of a real symmetric polynomial matrix.
The constant matrix must have distinct eigenvalues for a unique Taylor expansion.
Each eigenpair is lifted degree by degree with the constant bordered Jacobian.
"""
function eigh(A::AbstractMatrix{<:DA})
    lock(engine_lock) do
        foreach(valid, A)
        n, m = size(A)
        n == m || throw(DimensionMismatch("Expected a square matrix"))
        all(norm(A[i,j]-A[j,i], 0) == 0 for i in 1:n, j in 1:n) || throw(ArgumentError("Expected a symmetric matrix"))
        if all(iszero(A[i,j]) for i in 1:n for j in 1:n if i != j)
            permutation = sortperm(cons.(LinearAlgebra.diag(A)))
            return copy.(LinearAlgebra.diag(A)[permutation]), DA.(Matrix{Float64}(I, n, n)[:,permutation])
        end
        F = LinearAlgebra.eigen(LinearAlgebra.Symmetric(cons.(A)))
        values, vectors = DA.(F.values), DA.(F.vectors)
        isconstant = all(norm(trim(a, 1), 0) == 0 for a in A)
        isconstant && return values, vectors
        scale = max(1.0, maximum(abs, F.values; init=0.0))
        minimum((abs(F.values[i]-F.values[j]) for i in 1:n for j in i+1:n); init=Inf) > 100eps(Float64)*scale ||
            throw(ArgumentError("Repeated or numerically indistinguishable eigenvalues do not define unique Taylor eigenvectors"))
        saved = getTO()
        try
            for k in 1:n
                v0, lambda0 = F.vectors[:,k], F.values[k]
                J = [cons.(A)-lambda0*I -v0; permutedims(v0) 0.0]
                Ji = inv(J)
                v, lambda = copy(vectors[:,k]), values[k]
                for degree in 1:saved
                    setTO(degree)
                    residual = vcat(A*v - lambda*v, (sum(abs2, v)-1)/2)
                    correction = Ji*residual
                    v = v-correction[1:n]
                    lambda = lambda-correction[end]
                end
                vectors[:,k] = v
                values[k] = lambda
            end
        finally
            setTO(saved)
        end
        values, vectors
    end
end
function eigh(A::AbstractMatrix{<:Real})
    F = LinearAlgebra.eigen(LinearAlgebra.Symmetric(Matrix(A)))
    F.values, F.vectors
end
