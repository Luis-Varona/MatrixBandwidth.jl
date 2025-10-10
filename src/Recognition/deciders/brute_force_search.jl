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

# Supertype Hierarchy
`BruteForceSearch` <: [`AbstractDecider`](@ref) <: [`AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix, this brute-force algorithm runs in ``O(n! ⋅ n²)`` time:
- Up to ``n!/2`` permutations may be checked (except when ``n = 1``, in which case
    ``1! = 1`` permutation is checked). This is, clearly, ``O(n!)``.
- For each permutation, the [`bandwidth`](@ref) function is called on `view(A, perm, perm)`,
    which takes ``O(n²)`` time.
- Therefore, the overall time complexity is ``O(n! ⋅ n²)``.

# Examples
In many cases, the algorithm iterates over all (if ``k`` is smaller than the true minimu
bandwidth) or almost all (if ``k`` is equally to or only slightly larger than the true
minimum) possible permutations—in these cases, it is infeasible to go above ``9×9`` or
``10×10`` without incurring multiple-hour runtimes. (Even when ``k`` is considerably larger
than the true minimum, it is unlikely that a bandwidth-``k`` ordering will be found in a
reasonable time frame.) Nevertheless, we see that it is quite effective for, say, ``8×8``:
```@repl
using Random, SparseArrays
Random.seed!(314159);
(n, p) = (8, 0.5);
A = sprand(n, n, p);
A = A + A' # Ensure structural symmetry;
(k_false, k_true) = (3, 5);
has_bandwidth_k_ordering(A, k_false, Recognition.BruteForceSearch())
has_bandwidth_k_ordering(A, k_true, Recognition.BruteForceSearch())
```

# Notes
Brute force is by far the slowest approach to matrix bandwidth minimization and should only
be used in very niche cases (like verifying the correctness of other algorithms in unit
tests). For ``10×10`` matrices, the algorithm already takes several minutes to run for
difficult values of ``k`` (namely, values below or only slightly above the true minimum) and
allocates several gigabytes of memory. Given the ``O(n! ⋅ n²)`` time complexity, checking
"bandwidth ≤ k" would take over an hour for many ``11×11`` matrices.

This holds true even when ``k`` is considerably larger than the true minimum bandwidth—as
long as it remains below the bandwidth induced by the original ordering, it is unlikely that
a bandwidth-``k`` ordering will be found early simply by random chance. Additionally, time
complexity will remain on the order of ``n! ⋅ n²`` in the average case.

See also [`MatrixBandwidth.Minimization.Exact.BruteForceSearch`](@ref) for the minimization
variant of this algorithm (which simply never breaks early, instead iterating over all
permutations up to reversal to ensure that the minimum bandwidth is found).
"""
struct BruteForceSearch <: AbstractDecider end

push!(MatrixBandwidth.ALGORITHMS[:Recognition], BruteForceSearch)

Base.summary(::BruteForceSearch) = "Brute-force search"

MatrixBandwidth._requires_structural_symmetry(::BruteForceSearch) = false

#= We take advantage of the laziness of `permutations` and `Iterators.filter` to avoid
iterating over all orderings if a valid one is found early. =#
function _has_bandwidth_k_ordering_impl(
    A::AbstractMatrix{Bool}, k::Integer, ::BruteForceSearch
)
    #= `i₁, i₂, … iₙ` induces the same bandwidth as `iₙ, iₙ₋₁, … i₁`, so without loss of
    generality, we restrict our search to orderings such that `i₁ ≤ iₙ` (with equality
    checked just in case `n = 1`). =#
    orderings_up_to_reversal = Iterators.filter(
        perm -> perm[1] <= perm[end], permutations(axes(A, 1))
    )
    valid_orderings = Iterators.filter(
        perm -> bandwidth(view(A, perm, perm)) <= k, orderings_up_to_reversal
    )
    res = iterate(valid_orderings)

    if isnothing(res)
        ordering = nothing
    else
        ordering = res[1]
    end

    return ordering
end
