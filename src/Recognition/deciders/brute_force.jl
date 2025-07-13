# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    BruteForceSearch <: AbstractDecider <: AbstractAlgorithm

The simplest method for determining, given some fixed ``k ∈ ℕ``, whether a matrix has
bandwidth at most ``k`` up to symmetric permutation is to iterate over all orderings and
compute the bandwidth induced by each.

Since ``i₁, i₂, … iₙ`` induces the same bandwidth as ``iₙ, iₙ₋₁, … i₁``, we restrict our
search to orderings such that ``i₁ ≤ iₙ`` (with equality checked just in case ``n = 1``).

If a bandwidth-``k`` ordering is found, the algorithm breaks early instead of continuing
to check subsequent permutations.

# Performance
Given an ``n×n`` input matrix ``A``, this brute-force algorithm runs in ``O(n! ⋅ n²)`` time:
- Up to ``n!/2`` permutations may be checked (except when ``n = 1``, in which case
    ``1! = 1`` permutation is checked). This is, clearly, ``O(n!)``.
- For each permutation, the [`bandwidth`](@ref) function is called on ``A[perm, perm]``,
    which takes ``O(n²)`` time.
- Therefore, the overall time complexity is ``O(n! ⋅ n²)``.

# Examples
[TODO: Write here]

# Notes
Brute force is by far the slowest approach to matrix bandwidth recognition and should only
be used in very niche cases like writing unit tests for other non-naïve algorithms.
"""
struct BruteForceSearch <: AbstractDecider end

Base.summary(::BruteForceSearch) = "Brute-force search"

_requires_symmetry(::BruteForceSearch) = false

#= We take advantage of the laziness of `permutations` and `Iterators.filter` to avoid
iterating over all orderings if a valid one is found early. =#
function _bool_bandwidth_k_ordering(A::AbstractMatrix{Bool}, k::Int, ::BruteForceSearch)
    #= `i₁, i₂, … iₙ` induces the same bandwidth as `iₙ, iₙ₋₁, … i₁`, so without loss of
    generality, we restrict our search to orderings such that `i₁ ≤ iₙ` (with equality
    checked just in case `n = 1`). =#
    orderings_up_to_reversal = Iterators.filter(
        perm -> perm[1] <= perm[end], permutations(axes(A, 1))
    )
    valid_orderings = Iterators.filter(
        perm -> bandwidth(A[perm, perm]) <= k, orderings_up_to_reversal
    )
    res = iterate(valid_orderings)

    if isnothing(res)
        ordering = nothing
    else
        ordering = res[1]
    end

    return ordering
end
