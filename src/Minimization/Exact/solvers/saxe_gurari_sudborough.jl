# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SaxeGurariSudborough <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

The *Saxe–Gurari–Sudborough minimization algorithm* is an exact method for minimizing the
bandwidth of a structurally symmetric matrix ``A``. For a fixed ``k ∈ ℕ``, the algorithm
invokes a subroutine that determines whether ``A`` has bandwidth at most ``k`` up to
symmetric permutation. This subroutine employs dynamic programming to search over
equivalence classes of partial orderings, where two partial orderings of length ``l`` are
equivalent if they share the same *active region*. (The active region of a partial ordering
is defined as the sequence of the last ``min(k, l)`` vertices in the ordering taken together
with all *dangling edges*—edges with one endpoint in the ordering and the other endpoint not
yet in the ordering.) It extends these partial layouts one vertex at a time in a
breadth-first manner, pruning implausible classes that violate bandwidth-``k`` constraints
such as degree bounds on active vertices and excessive numbers of dangling edges [GS84].
This search is repeated with incrementing values of ``k`` until a bandwidth-``k`` ordering
is found, with ``k`` initialized to some lower bound on the minimum bandwidth of ``A`` up to
symmetric permutation.

Specifically, this implementation of the Saxe–Gurari–Sudborough algorithm uses the
``min(α(A), γ(A))`` lower bound from [CS05, pp. 359--60] as the initial value of ``k``.
(Further implementation details can be found in the source code for
[`bandwidth_lower_bound`](@ref).)

As noted above, the Saxe–Gurari–Sudborough algorithm requires structurally symmetric input
(that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for
``1 ≤ i, j ≤ n``).

# Supertype Hierarchy
`SaxeGurariSudborough` <: [`ExactSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix, the Saxe–Gurari–Sudborough algorithm runs in ``O(nⁿ⁻¹)``
time:
- For each underlying "bandwidth ≤ k" check, we call the Saxe–Gurari–Sudborough recognition
  algorithm, which runs in ``O(nᵏ)`` time [GS84, p. 531]. (This is an improvement upon the
  original ``O(nᵏ⁺¹)`` Saxe algorithm [Sax80, p. 363].)
- The difference between the maximum possible bandwidth (``n - 1``) and our initial lower
    bound grows linearly in ``n``, so we run the underlying ``O(nᵏ)`` recognition algorithm
    ``O(n)`` times.
- Therefore, the overall time complexity is ``∑ₖ₌₀ⁿ⁻¹ nᵏ = O(nⁿ⁻¹)``.

Whereas most exact bandwidth minimization algorithms are technically factorial-time (with
respect to ``n``) in the worst case but practically always approximate exponential time
complexity in real life, the ``O(nⁿ⁻¹)`` upper bound on Saxe–Gurari–Sudborough is typically
a good representation of actual performance in most cases. Indeed, these other types of
algorithms tend to outperform Saxe–Gurari–Sudborough for larger ``n``, given that their
aggressive pruning strategies keep their effective search space very small in practice and
``O(2ⁿ)`` ⊂ ``O(nⁿ⁻¹)``.

# Examples
We verify the optimality of the ordering found by Saxe–Gurari–Sudborough for a random
``9×9`` matrix via a brute-force search over all possible permutations up to reversal:
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(52452);

julia> (n, p) = (9, 0.5);

julia> A = sprand(n, n, p);

julia> A = A + A' # Ensure structural symmetry;

julia> res_bf = minimize_bandwidth(A, Minimization.BruteForceSearch())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Brute-force search
 * Approach: exact
 * Minimum Bandwidth: 5
 * Original Bandwidth: 8
 * Matrix Size: 9×9

julia> res_sgs = minimize_bandwidth(A, Minimization.SaxeGurariSudborough())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Saxe–Gurari–Sudborough
 * Approach: exact
 * Minimum Bandwidth: 5
 * Original Bandwidth: 8
 * Matrix Size: 9×9
```

We now generate (and shuffle) a random ``25×25`` matrix with minimum bandwidth ``5`` using
[`MatrixBandwidth.random_banded_matrix`](@ref). Saxe–Gurari–Sudborough then finds a
bandwidth-``4`` ordering, which is (we claim) optimal up to symmetric permutation. (In some
cases, `random_banded_matrix(n, k)` generates matrices with minimum bandwidth `< k`—this
appears to be one such case. Although we do not explicitly verify exact optimality—which
*is* guaranteed by the original paper [GS84]—here via brute-force search, this example
demonstrates that Saxe–Gurari–Sudborough at the very least finds a quite good ordering.)
```jldoctest
julia> using Random

julia> Random.seed!(937497);

julia> (n, k, p) = (25, 5, 0.25);

julia> A = MatrixBandwidth.random_banded_matrix(n, k; p=p);

julia> perm = randperm(n);

julia> A_shuffled = A[perm, perm];

julia> bandwidth(A)
5

julia> bandwidth(A_shuffled) # Much larger after shuffling
19

julia> minimize_bandwidth(A_shuffled, Minimization.SaxeGurariSudborough())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Saxe–Gurari–Sudborough
 * Approach: exact
 * Minimum Bandwidth: 4
 * Original Bandwidth: 19
 * Matrix Size: 25×25
```

# Notes
The Saxe–Gurari–Sudborough algorithm was originally a bandwidth recognition algorithm, not a
minimization one—as previously mentioned, we repurpose it here by repeatedly invoking the
original procedure for incrementing values of ``k`` until a valid ordering is found. The
general family of recognition algorithms to which Saxe–Gurari–Sudborough belongs was
conceived as a response to a question posed by [GGJK78, p. 494]: is the "bandwidth ≤ k?"
problem NP-complete for arbitrary ``k``? [Sax80] answered this question in the negative by
providing a ``O(nᵏ⁺¹)`` algorithm, constructively proving that the problem is class P.
Later, [GS84] improved upon this algorithm by reducing time complexity to ``O(nᵏ)``. Whereas
the original Saxe algorithm considers extensions of partial orderings with any remaining
unplaced vertex (of which there are ``O(n)`` at any point in the breadth-first search), the
Gurari–Sudborough refinement only considers extensions with vertices reachable by paths
beginning with a dangling edge that never again traverse a dangling edge [GS84, pp.
535–36].

# References
- [GGJK78](@cite): M. R. Garey, R. L. Graham, D. S. Johnson and D. E. Knuth. *Complexity
    Results for Bandwidth Minimization*. SIAM Journal on Applied Mathematics **34**, 477–95
    (1978). https://doi.org/10.1137/0134037.
- [GS84](@cite): E. M. Gurari and I. H. Sudborough. *Improved dynamic programming algorithms
    for bandwidth minimization and the MinCut Linear Arrangement problem*. Journal of
    Algorithms **5**, 531–46 (1984). https://doi.org/10.1016/0196-6774(84)90006-3.
- [Sax80](@cite): J. B. Saxe. *Dynamic-Programming Algorithms for Recognizing
    Small-Bandwidth Graphs in Polynomial Time*. SIAM Journal on Algebraic and Discrete
    Methods **1**, 363–69 (1980). https://doi.org/10.1137/0601042.
"""
struct SaxeGurariSudborough <: ExactSolver end

push!(MatrixBandwidth.ALGORITHMS[:Minimization][:Exact], SaxeGurariSudborough)

Base.summary(::SaxeGurariSudborough) = "Saxe–Gurari–Sudborough"

MatrixBandwidth._requires_structural_symmetry(::SaxeGurariSudborough) = true

function Minimization._minimize_bandwidth_impl(
    A::AbstractMatrix{Bool}, ::SaxeGurariSudborough
)
    components = connected_components(A)
    #= Heuristically, it is likelier that smaller components have lower bandwidth, so we
    process them first to keep `k` low for as long as possible (since the complexity of each
    recognition subroutine is exponential in `k`). =#
    sort!(components; by=length)

    ordering = Vector{Int}(undef, size(A, 1))
    k = 1
    num_placed = 0

    for component in components
        submatrix = view(A, component, component)
        k = max(k, bandwidth_lower_bound(submatrix))
        component_ordering = Recognition._sgs_connected_ordering(submatrix, k)

        while isnothing(component_ordering)
            component_ordering = Recognition._sgs_connected_ordering(submatrix, k += 1)
        end

        ordering[(num_placed + 1):(num_placed += length(component))] .= component[component_ordering]
    end

    return ordering
end
