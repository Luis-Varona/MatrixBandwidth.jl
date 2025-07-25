# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    BruteForceSearch <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

The simplest exact method for minimizing the bandwidth of a matrix is to iterate over all
possible symmetric permutations and compare the bandwidths they induce.

Since ``i₁, i₂, … iₙ`` induces the same bandwidth as ``iₙ, iₙ₋₁, … i₁``, we restrict our
search to orderings such that ``i₁ ≤ iₙ`` (with equality checked just in case ``n = 1``).

# Performance
Given an ``n×n`` input matrix ``A``, this brute-force algorithm runs in ``O(n! ⋅ n²)`` time:
- Precisely ``n!/2`` permutations are checked (except when ``n = 1``, in which case
    ``1! = 1`` permutation is checked). This is, clearly, ``O(n!)``.
- For each permutation, the [`bandwidth`](@ref) function is called on `view(A, perm, perm)`,
    which takes ``O(n²)`` time.
- Therefore, the overall time complexity is ``O(n! ⋅ n²)``.

Indeed, due to the need to exhaustively check all permutations, this is close to a lower
bound as well on the the algorithm's time complexity. (The only reason we cannot claim to
have a precise value for the big-``Θ`` complexity is that the [`bandwidth`](@ref) function
is not *exactly* ``Θ(n²)``, although it is close.)

# Examples
The algorithm always iterates over all possible permutations, so it is infeasible to go
above ``9×9`` or ``10×10`` without incurring multiple-hour runtimes. Nevertheless, we see
that it is quite effective for, say, ``8×8``:
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(628318);

julia> (n, p) = (8, 0.2);

julia> A = sprand(Bool, n, n, p);

julia> A = A + A' # Ensure structural symmetry;

julia> minimize_bandwidth(A, Minimization.BruteForceSearch())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Brute-force search
 * Approach: exact
 * Minimum Bandwidth: 3
 * Original Bandwidth: 6
 * Matrix Size: 8×8
```

# Notes
Brute force is by far the slowest approach to matrix bandwidth minimization and should only
be used in very niche cases (like verifying the correctness of other algorithms in unit
tests). For ``10×10`` matrices, the algorithm already takes several minutes to run (between
``2`` to ``5`` on most commercial machines) and allocates over ``4`` gigabytes of memory.
Given the ``O(n! ⋅ n²)`` time complexity, minimizing the bandwidth of any ``11×11`` matrix
would take over an hour.

"""
struct BruteForceSearch <: ExactSolver end

Base.summary(::BruteForceSearch) = "Brute-force search"

_requires_symmetry(::BruteForceSearch) = false

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::BruteForceSearch)
    #= `i₁, i₂, … iₙ` induces the same bandwidth as `iₙ, iₙ₋₁, … i₁`, so without loss of
    generality, we restrict our search to orderings such that `i₁ ≤ iₙ` (with equality
    checked just in case `n = 1`). =#
    orderings_up_to_reversal = Iterators.filter(
        perm -> perm[1] <= perm[end], permutations(axes(A, 1))
    )
    return argmin(perm -> bandwidth(view(A, perm, perm)), orderings_up_to_reversal)
end
