# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    has_bandwidth_k_ordering(A, k, decider=CapraraSalazarGonzalez()) -> BandRecogResult

Determine whether `A` has bandwidth at most `k` using the algorithm defined by `decider`.

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

This function [TODO: Write here]

# Arguments
- `A::AbstractMatrix{T}`: the (square) matrix whose bandwidth is to be tested.
- `k::Int`: the threshold bandwidth against which to test.
- `decider::AbstractDecider`: the matrix bandwidth recognition algorithm to use; defaults to
    [`CapraraSalazarGonzalez`](@ref). (See the [`Recognition`](@ref) module documentation
    for a full list of supported deciders.)

# Returns
- `::BandRecogResult`: TODO: Write here

# Examples
[TODO: Add here once more deciders are implemented]

# Notes
Some texts define matrix bandwidth to be the minimum non-negative integer ``k`` such that
``A[i, j] = 0`` whenever ``|i - j| ≥ k`` instead, particularly in more mathematically-minded
communities. Effectively, this definition treats diagonal matrices as bandwidth ``1``,
tridiagonal matrices as bandwidth ``2``, and so on. Our definition, on the other hand, is
more common in computer science contexts, treating diagonal matrices as bandwidth ``0`` and
tridiagonal matrices as bandwidth ``1``. (Both definitions, however, agree that the
bandwidth of an empty matrix is simply ``0``.)
"""
function has_bandwidth_k_ordering(
    A::AbstractMatrix{T}, k::Int, decider::AbstractDecider=DEFAULT_DECIDER
) where {T<:Number}
    # TODO: Implement
end
