# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    bandwidth(A) -> Int

Compute the bandwidth of `A` without any permutations.

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

In contrast to [`minimize_bandwidth`](@ref), this function does not attempt to minimize the
bandwidth of `A` by permuting its rows and columns—it simply computes its bandwidth as is.

# Arguments
- `A::AbstractMatrix{T}`: the (square) matrix whose bandwidth is to be computed.

# Returns
- `::Int`: the bandwidth of `A`.

# Examples
`bandwidth` correctly identifies the bandwidth of a pentadiagonal matrix as ``2`` and does
not attempt to find a minimizing permutation upon shuffling of its rows and columns:
```jldoctest
julia> using Random

julia> Random.seed!(242622);

julia> (n, k) = (8, 2);

julia> perm = randperm(n);

julia> A = (!iszero).(random_banded_matrix(8, 2))
8×8 BitMatrix:
 1  0  1  0  0  0  0  0
 1  0  1  1  0  0  0  0
 1  1  1  1  1  0  0  0
 0  1  0  1  1  1  0  0
 0  0  0  1  1  1  1  0
 0  0  0  1  1  1  1  1
 0  0  0  0  1  1  1  1
 0  0  0  0  0  1  0  1

julia> bandwidth(A)
2

julia> A_shuffled = A[perm, perm]
8×8 BitMatrix:
 1  1  1  0  0  0  1  1
 1  1  1  0  0  0  0  1
 1  1  1  0  0  0  1  0
 0  0  1  1  1  1  1  0
 0  0  0  1  1  0  0  0
 0  0  0  1  1  0  1  0
 1  0  1  0  0  1  1  0
 1  0  0  0  0  0  0  1

julia> bandwidth(A_shuffled)
7
```

# Notes
Some texts define matrix bandwidth to be the minimum non-negative integer ``k`` such that
``A[i, j] = 0`` whenever ``|i - j| ≥ k`` instead, particularly in more mathematically-minded
communities. Effectively, this definition treats diagonal matrices as bandwidth ``1``,
tridiagonal matrices as bandwidth ``2``, and so on. Our definition, on the other hand, is
more common in computer science contexts, treating diagonal matrices as bandwidth ``0`` and
tridiagonal matrices as bandwidth ``1``. (Both definitions, however, agree that the
bandwidth of an empty matrix is simply ``0``.)
"""
function bandwidth(A::AbstractMatrix{T}) where {T<:Number}
    _assert_matrix_is_square(A)

    if T === Bool
        A_bool = A
    else
        A_bool = (!iszero).(A)
    end

    indices = findall(A_bool)

    if isempty(indices)
        band = 0
    else
        band = maximum(abs(index[1] - index[2]) for index in indices)
    end

    return band
end
