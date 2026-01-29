# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CapraraSalazarGonzalez <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

The *Caprara–Salazar-González minimization algorithm* is an exact method for minimizing the
bandwidth of a structurally symmetric matrix ``A``. For a fixed ``k ∈ ℕ``, the algorithm
performs a depth-first search of all partial orderings of the rows and columns of ``A``,
adding indices one at a time. Partial orderings are pruned not only by ensuring that
adjacent pairs of currently placed indices are within ``k`` of each other but also by
employing a branch-and-bound framework with lower bounds on bandwidth compatibility computed
via integer linear programming relaxations. This search is repeated with incrementing values
of ``k`` until a bandwidth-``k`` ordering is found [CS05], with ``k`` initialized to some
lower bound on the minimum bandwidth of ``A`` up to symmetric permutation.

Specifically, this implementation of the Caprara–Salazar-González algorithm uses the
``min(α(A), γ(A))`` lower bound from the original paper [CS05, pp. 359--60] as the initial
value of ``k``. (Further implementation details can be found in the source code for
[`bandwidth_lower_bound`](@ref).)

As noted above, the Caprara–Salazar-González algorithm requires structurally symmetric input
(that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for
``1 ≤ i, j ≤ n``).

# Supertype Hierarchy
`CapraraSalazarGonzalez` <: [`ExactSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix with ``m`` nonzero entries, the Caprara–Salazar-González
algorithm runs in ``O(n! ⋅ mn²)`` time in the worst case:

- For each underlying "bandwidth ≤ ``k``" check, we perform a depth-first search of
    ``O(n!)`` partial orderings.
- At each search node, we compute first and last feasible positions for all free nodes in
    ``O(mn)`` time using the closed-form formulas from Propositions 11 and 14 of [CS05].
    Thus, the overall time complexity for each value of ``k`` is ``O(n! ⋅ mn)``.
- The difference between the maximum possible bandwidth (``n - 1``) and our initial lower
    bound grows linearly in ``n``, so we run the underlying ``O(n! ⋅ mn)`` recognition
    algorithm ``O(n)`` times.
- Therefore, the overall time complexity is ``O(n! ⋅ mn ⋅ n) = O(n! ⋅ mn²)``.

Of course, this is all but an upper bound on the time complexity of
Caprara–Salazar-González, achieved only in the most pathological of cases. In practice,
efficient pruning techniques and compatibility checks—along with [CS05, pp. 359--60]'s
relatively tight initial lower bound on the minimum bandwidth—result in approximately
exponential growth in time complexity with respect to ``n``.

# Examples
We verify the optimality of the ordering found by Caprara–Salazar-González for a random
``8×8`` matrix via a brute-force search over all possible permutations up to reversal:
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(5883);

julia> (n, p) = (8, 0.25);

julia> A = sprand(n, n, p);

julia> A = A + A' # Ensure structural symmetry;

julia> res_bf = minimize_bandwidth(A, Minimization.BruteForceSearch())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Brute-force search
 * Approach: exact
 * Minimum Bandwidth: 4
 * Original Bandwidth: 7
 * Matrix Size: 8×8

julia> res_csg = minimize_bandwidth(A, Minimization.CapraraSalazarGonzalez())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Caprara–Salazar-González
 * Approach: exact
 * Minimum Bandwidth: 4
 * Original Bandwidth: 7
 * Matrix Size: 8×8
```

# Notes
For readers of the original paper, what we call the Caprara–Salazar-González algorithm here
is designated the `LAYOUT_LEFT_TO_RIGHT` algorithm in [CS05]. The paper also describes a
`LAYOUT_BOTH_WAYS` algorithm that performs a bidirectional search by adding indices to both
the left and right ends of the current partial ordering. However, this version is
considerably more complex to implement, and we ran into problems enforcing ILP constraints
on node pairs added to opposite ends of the ordering. In any case, computational results
demonstrate that neither `LAYOUT_LEFT_TO_RIGHT` nor `LAYOUT_BOTH_WAYS` is consistently
faster, and the paper states that there is no known heuristic for determining which version
will be more performant for a given input [CS05, pp. 368--69]. Therefore, we opt to
implement only `LAYOUT_LEFT_TO_RIGHT` as a matter of practicality, although future
developers may wish to extend the interface with `LAYOUT_BOTH_WAYS` as well.

# References
- [CS05](@cite): A. Caprara and J.-J. Salazar-González. *Laying Out Sparse Graphs with
    Provably Minimum Bandwidth*. INFORMS Journal on Computing **17**, 356–73 (2005).
    https://doi.org/10.1287/ijoc.1040.0083.
"""
struct CapraraSalazarGonzalez <: ExactSolver end

push!(MatrixBandwidth.ALGORITHMS[:Minimization][:Exact], CapraraSalazarGonzalez)

Base.summary(::CapraraSalazarGonzalez) = "Caprara–Salazar-González"

MatrixBandwidth._requires_structural_symmetry(::CapraraSalazarGonzalez) = true

function Minimization._minimize_bandwidth_impl(
    A::AbstractMatrix{Bool}, ::CapraraSalazarGonzalez
)
    n = size(A, 1)

    ordering_buf = Vector{Int}(undef, n)
    k = bandwidth_lower_bound(A)
    dist_matrix = floyd_warshall_shortest_paths(A)
    fixed = Int[]
    unselected = Set(1:n)

    ordering = nothing

    while isnothing(ordering)
        ordering = Recognition._csg_layout_left_to_right!(
            ordering_buf, A, k, dist_matrix, fixed, unselected
        )
        k += 1
    end

    return ordering
end
