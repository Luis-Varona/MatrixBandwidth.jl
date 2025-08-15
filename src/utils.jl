# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    random_banded_matrix(n, k; p=0.5, rng=default_rng()) -> Matrix{Float64}

Generate a random `n×n` structurally symmetric `k`-banded matrix with band density `≈ p`.

By definition of structural symmetry, the ``(i, j)``-th entry of the matrix is nonzero if
and only if the ``(j, i)``-th entry is nonzero as well. All entries from this matrix are
from the interval `[0, 1]`. Entries up to the `k`ᵗʰ superdiagonal and down to the `k`ᵗʰ
subdiagonal are nonzero with probability `p`.

It is also guaranteed that each of these bands (besides the main diagonal) has at least one
nonzero entry (even when `p` is very small), thus ensuring that the matrix has bandwidth
precisely `k` before any reordering. (There may, however, still exist a symmetric
permutation inducing a minimum bandwidth less than `k`, especially for small values of `p`.)

# Arguments
- `n::Integer`: the order of the matrix to generate. Must be positive.
- `k::Integer`: the desired matrix bandwidth. Must satisfy `0 ≤ k < n`.

# Keyword Arguments
- `p::Real=0.5`: the band density. Must satisfy `0 < p ≤ 1`. Defaults to `0.5`.
- `rng::AbstractRNG=Random.default_rng()`: the random number generator to use. Defaults to
    `Random.default_rng()`.

# Returns
- `::Matrix{Float64}`: a random `n×n` matrix with bandwidth exactly `k` and sparse bands
    with density `p`.

# Examples
Generate a ``6×6`` matrix with bandwidth ``1`` and the maximum number of nonzero entries:
```jldoctest
julia> using Random

julia> A = random_banded_matrix(6, 1; p=1, rng=MersenneTwister(1228))
6×6 Matrix{Float64}:
 0.310239  0.346413  0.0       0.0        0.0       0.0
 0.509981  0.917073  0.390771  0.0        0.0       0.0
 0.0       0.760045  0.808396  0.0195686  0.0       0.0
 0.0       0.0       0.222338  0.853164   0.806888  0.0
 0.0       0.0       0.0       0.421603   0.132165  0.805813
 0.0       0.0       0.0       0.0        0.305339  0.0799183

julia> bandwidth(A)
1
```

Generate a ``7×7`` matrix with bandwidth ``3`` and band density ``0.3``:
```jldoctest
julia> using Random

julia> A = random_banded_matrix(7, 3; p=0.3, rng=MersenneTwister(0402))
7×7 Matrix{Float64}:
 0.0       0.132699  0.0       0.0       0.0  0.0       0.0
 0.869352  0.0       0.324319  0.926496  0.0  0.0       0.0
 0.0       0.891878  0.0       0.658102  0.0  0.0       0.0
 0.0       0.88859   0.399559  0.0       0.0  0.284285  0.703377
 0.0       0.0       0.0       0.0       0.0  0.0       0.0
 0.0       0.0       0.0       0.489594  0.0  0.0       0.393573
 0.0       0.0       0.0       0.412412  0.0  0.47063   0.0

julia> bandwidth(A)
3
```

Generate an ``8×8`` diagonal (bandwidth ``0``) matrix with default band density (``0.5``):
```jldoctest
julia> using Random

julia> A = random_banded_matrix(8, 0; rng=MersenneTwister(0102))
8×8 Matrix{Float64}:
 0.0  0.0        0.0       0.0       0.0  0.0      0.0  0.0
 0.0  0.0762399  0.0       0.0       0.0  0.0      0.0  0.0
 0.0  0.0        0.373113  0.0       0.0  0.0      0.0  0.0
 0.0  0.0        0.0       0.726309  0.0  0.0      0.0  0.0
 0.0  0.0        0.0       0.0       0.0  0.0      0.0  0.0
 0.0  0.0        0.0       0.0       0.0  0.41974  0.0  0.0
 0.0  0.0        0.0       0.0       0.0  0.0      0.0  0.0
 0.0  0.0        0.0       0.0       0.0  0.0      0.0  0.293132

julia> bandwidth(A)
0
```

# Notes
Users of the [`MatrixBandwidth`](@ref) package may find this function useful when generating
random test data for whatever frameworks, algorithms, etc. they are implementing.
"""
function random_banded_matrix(
    n::Integer, k::Integer; p::Real=0.5, rng::AbstractRNG=Random.default_rng()
)
    if n <= 0
        throw(ArgumentError("Matrix order must be positive, got $n"))
    end

    if k < 0
        throw(ArgumentError("Matrix bandwidth must be non-negative, got $k"))
    end

    if n <= k
        throw(ArgumentError("Matrix order must be greater than bandwidth, got $n and $k"))
    end

    if !(0 < p <= 1)
        throw(ArgumentError("Band density must be in (0, 1], got $p"))
    end

    A = zeros(Float64, n, n)

    #= The main diagonal is handled separately to halve the number of operations and avoid
    the need to ensure at least one nonzero entry (as this would not affect bandwidth). =#
    for i in 1:n
        if rand(rng) < p
            A[i, i] = rand(rng)
        end
    end

    #= Conditionally add entries to the superdiagonals and subdiagonals (up to the `k`ᵗʰ)
    based on the band density `p`. =#
    for d in 1:k
        feasible_indices = 1:(n - d)

        for i in feasible_indices
            if rand(rng) < p
                A[i, i + d] = rand(rng)
                A[i + d, i] = rand(rng)
            end
        end

        #= Ensure that the `d`ᵗʰ superdiagonal and subdiagonal each have at least one
        nonzero entry for `1 ≤ d ≤ k`, making a `bandwidth < k` symmetric permutation less
        likely. (Given structural symmetry, we need only check for empty superdiagonals.) =#
        if all(iszero, Iterators.map(i -> A[i, i + d], feasible_indices))
            i = rand(rng, feasible_indices)
            A[i, i + d] = rand(rng)
            A[i + d, i] = rand(rng)
        end
    end

    return A
end

# Find the indices of all connected components in an adjacency matrix
function _connected_components(A::AbstractMatrix{Bool})
    n = size(A, 1)
    visited = falses(n)
    queue = Queue{Int}()
    components = Vector{Int}[]

    for i in 1:n
        if !visited[i]
            visited[i] = true
            enqueue!(queue, i)
            component = Int[]

            while !isempty(queue)
                u = dequeue!(queue)
                push!(component, u)

                for v in findall(view(A, :, u))
                    if !visited[v]
                        visited[v] = true
                        enqueue!(queue, v)
                    end
                end
            end

            push!(components, component)
        end
    end

    return components
end

# Identify the highest supertype of `subtype` that is a subtype of `abstracttype`
function _find_direct_subtype(abstracttype::Type, subtype::Type)
    if !isabstracttype(abstracttype)
        throw(ArgumentError("Expected an abstract type, got $abstracttype"))
    end

    if !(subtype <: abstracttype)
        throw(ArgumentError("Expected a subtype of $abstracttype, got $subtype"))
    end

    parent = supertype(subtype)

    if parent === abstracttype
        child = subtype
    else
        child = _find_direct_subtype(abstracttype, parent)
    end

    return child
end

#= Check whether `A[i, j]` is nonzero if and only if `A[j, i]` is nonzero for all `i` and
`j`. Instead of transposing `A` and allocating a new matrix, it suffices to iterate over all
opposing pairs of off-diagonal entries. (`A` is assumed to be square.) =#
function _is_structurally_symmetric(A::AbstractMatrix{<:Number})
    n = size(A, 1)
    return all((A[i, j] != 0) == (A[j, i] != 0) for i in 1:(n - 1) for j in (i + 1):n)
end

#= Convert `A` to a  to improve performance via cache optimizations, bitwise
operations, etc. Any symmetric permutation preserves diagonal entries, so the diagonal of
the output matrix is set to `false` for consistency with any algorithms that require an
adjacency matrix structure. =#
function _offdiag_nonzero_support(A::AbstractMatrix{T}) where {T<:Number}
    if T === Bool
        A_bool = BitMatrix(A) # Copy to avoid shared mutability
    else
        A_bool = (!iszero).(A)
    end

    foreach(i -> A_bool[i, i] = false, axes(A_bool, 1))

    return A_bool
end
