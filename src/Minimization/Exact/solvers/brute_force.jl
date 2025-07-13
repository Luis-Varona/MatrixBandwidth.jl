# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    BruteForce <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

The simplest exact method for minimizing the bandwidth of a matrix is to iterate over all
possible symmetric permutations and compare the bandwidths they induce.

Since ``i₁, i₂, … iₙ`` induces the same bandwidth as ``iₙ, iₙ₋₁, … i₁``, we restrict our
search to orderings such that ``i₁ ≤ iₙ`` (with equality checked just in case ``n = 1``).

# Performance
Given an ``n×n`` input matrix ``A``, this brute-force algorithm runs in ``O(n! ⋅ n²)`` time:
- Precisely ``n!/2`` permutations are checked (except when ``n = 1``, in which case
    ``1! = 1`` permutation is checked). This is, clearly, ``O(n!)``.
- For each permutation, the [`bandwidth`](@ref) function is called on ``A[perm, perm]``,
    which takes ``O(n²)`` time.
- Therefore, the overall time complexity is ``O(n! ⋅ n²)``.

# Examples
[TODO: Write here]

# Notes
Brute force is by far the slowest approach to matrix bandwidth minimization and should only
be used in very niche cases like writing unit tests for other non-naïve algorithms.
"""
struct BruteForce <: ExactSolver end

Base.summary(::BruteForce) = "Brute-force search"

_requires_symmetry(::BruteForce) = false

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::BruteForce)
    #= `i₁, i₂, … iₙ` induces the same bandwidth as `iₙ, iₙ₋₁, … i₁`, so without loss of
    generality, we restrict our search to orderings such that `i₁ ≤ iₙ` (with equality
    checked just in case `n = 1`). =#
    orderings_up_to_reversal = Iterators.filter(
        perm -> perm[1] <= perm[end], permutations(axes(A, 1))
    )
    return argmin(perm -> bandwidth(A[perm, perm]), orderings_up_to_reversal)
end
