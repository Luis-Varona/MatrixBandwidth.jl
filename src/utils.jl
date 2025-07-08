# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    random_banded_matrix(n, k; p=0.75, rng=default_rng()) -> Matrix{Float64}

Generate a random `n×n` matrix with bandwidth exactly `k` and sparse bands with density `p`.

All entries from this matrix will be from the interval `[0, 1]`. Entries up to the `k`-th
superdiagonal and down to the `k`-th subdiagonal are nonzero with probability `p`, and each
band has at least one nonzero entry to ensure that the bandwidth is precisely `k`.

# Arguments
- `n::Int`: the order of the matrix to generate. Must be positive.
- `k::Int`: the desired matrix bandwidth. Must satisfy `0 ≤ k < n`.

# Keyword Arguments
- `p::Real=0.75`: the band density. Must satisfy `0 < p ≤ 1`.
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
 0.918835  0.816296   0.0       0.0        0.0       0.0
 0.182127  0.782844   0.616169  0.0        0.0       0.0
 0.0       0.0445171  0.916205  0.730272   0.0       0.0
 0.0       0.0        0.966811  0.414062   0.210912  0.0
 0.0       0.0        0.0       0.0150353  0.135984  0.558082
 0.0       0.0        0.0       0.0        0.428772  0.329567

julia> bandwidth(A)
1
```

Generate a ``7×7`` matrix with bandwidth ``3`` and band density `0.3`:
```jldoctest
julia> using Random

julia> A = random_banded_matrix(7, 3; p=0.3, rng=MersenneTwister(0402))
7×7 Matrix{Float64}:
 0.856072  0.720893  0.0       0.0       0.0       0.0        0.0
 0.0       0.0       0.0       0.646516  0.845229  0.0        0.0
 0.997473  0.773515  0.854375  0.926462  0.21636   0.0        0.0
 0.0       0.516052  0.220979  0.844818  0.0       0.0395003  0.568892
 0.0       0.402696  0.499802  0.0       0.304168  0.237423   0.0
 0.0       0.0       0.0       0.0       0.0       0.0        0.877917
 0.0       0.0       0.0       0.101071  0.0       0.0        0.829221

julia> bandwidth(A)
3
```

Generate an ``8×8`` diagonal (bandwidth ``0``) matrix with default band density (``0.75``):
```jldoctest
julia> using Random

julia> A = random_banded_matrix(8, 0; rng=MersenneTwister(0102))
8×8 Matrix{Float64}:
 0.781618  0.0      0.0       0.0  0.0        0.0      0.0       0.0
 0.0       0.56589  0.0       0.0  0.0        0.0      0.0       0.0
 0.0       0.0      0.966643  0.0  0.0        0.0      0.0       0.0
 0.0       0.0      0.0       0.0  0.0        0.0      0.0       0.0
 0.0       0.0      0.0       0.0  0.0412729  0.0      0.0       0.0
 0.0       0.0      0.0       0.0  0.0        0.70196  0.0       0.0
 0.0       0.0      0.0       0.0  0.0        0.0      0.494806  0.0
 0.0       0.0      0.0       0.0  0.0        0.0      0.0       0.227507

julia> bandwidth(A)
0
```

# Notes
Users of the [`MatrixBandwidth`](@ref) package may find this function useful when generating
random test data for whatever frameworks, algorithms, etc. they are implementing.
"""
function random_banded_matrix(
    n::Int, k::Int; p::Real=0.75, rng::AbstractRNG=Random.default_rng()
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

    # Ensure that the main diagonal has at least one nonzero entry
    diag_idx = rand(rng, 1:n)
    A[diag_idx, diag_idx] = rand(rng)

    #= Ensure that each superdiagonal and subdiagonal has at least one nonzero entry. (The
    `is_upper` flag indicates whether we are filling a superdiagonal or subdiagonal.) =#
    for offset in 1:k, is_upper in (true, false)
        rows = ((1 + offset):n) .- is_upper * offset
        cols = (1:(n - offset)) .+ is_upper * offset
        idx = rand(rng, 1:length(rows))
        A[rows[idx], cols[idx]] = rand(rng)
    end

    for i in 1:n, j in max(1, i - k):(min(n, i + k))
        if rand(rng) < p
            A[i, j] = rand(rng)
        end
    end

    return A
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
