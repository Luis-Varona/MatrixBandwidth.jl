# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    const RCM_INPUT :: BitMatrix

A sample input for testing the Cuthill–McKee and reverse Cuthill–McKee algorithms.

The original matrix bandwidth is 8, and the expected outputs are provided in
[`RCM_ANSWERS`](@ref). This example is taken from [LLS+01](@cite) (converted from C++'s
zero-based indexing to Julia's one-based indexing).
"""
const RCM_INPUT = let
    n = 10
    edges = [
        (1, 4),
        (1, 6),
        (2, 3),
        (2, 5),
        (2, 7),
        (2, 10),
        (3, 4),
        (3, 5),
        (4, 6),
        (4, 9),
        (5, 7),
        (6, 7),
        (6, 8),
        (7, 8),
    ]
    A = falses(n, n)
    foreach(e -> A[e...] = A[reverse(e)...] = true, edges)
    A
end

"""
    const RCM_ANSWERS :: Vector{Tuple{Int, Vector{Int}}}

Examples of outputs when reverse Cuthill–McKee is correctly applied to [`RCM_INPUT`](@ref).

The first element of each tuple is the bandwidth of the reordered matrix, and the second is
the new ordering of the rows and columns. This example is taken from [LLS+01](@cite)
(converted from C++'s zero-based indexing to Julia's one-based indexing).

Note that these are only some possible outputs, as different node selection heuristics may
yield different starting points and thus different orderings. Given matching starting nodes,
however, these are the expected results.

If testing Cuthill–McKee instead of reverse Cuthill–McKee, the orderings can simply be
reversed, although bandwidths may differ in the general case.
"""
const RCM_ANSWERS = [
    (4, [9, 4, 1, 10, 3, 6, 2, 5, 8, 7]),
    (4, [10, 2, 5, 7, 8, 3, 9, 6, 4, 1]),
    (4, [1, 9, 6, 8, 4, 7, 5, 3, 2, 10]),
]

"""
    random_banded_discon_matrix(n, k, p, num_ccs)

Generate a random disconnected `n×n` matrix with bandwidth `k`.

Specifically, this function generates `num_ccs` separate matrices with bandwidth precisely
`k` using [`random_banded_matrix`](@ref) then combines them into a single block diagonal
matrix. Treating the resulting matrix as the adjacency of a directed graph, this ensures
that there are no edges between any of the blocks, guaranteeing that the resultant graph is
disconnected with at least `num_ccs` connected components. (If any of the blocks themselves
are disconnected, the graph may have even more than `num_ccs` components.)

Since the Cuthill–McKee and reverse Cuthill–McKee algorithms are designed to work on
connected components (or their matrix representations), this acts as a useful stress test
when validating that our implementation extends to disconnected matrices.

# Arguments
- `n::Int`: the order of the matrix to generate. Must be positive.
- `k::Int`: the desired matrix bandwidth. Must satisfy `0 ≤ k < n`.
- `num_ccs::Int`: the minimum number of connected components. Must satisfy `1 ≤ num_ccs ≤ n`.
- `p::Real`: the band density. Must satisfy `0 < p ≤ 1`.

# Returns
- `::SparseMatrixCSC{Float64, Int}`: a random `n×n` matrix with bandwidth exactly `k` and at
    least `num_ccs` connected components.
"""
function random_banded_discon_matrix(n::Int, k::Int, num_ccs::Int, p::Real)
    if n <= 0
        throw(ArgumentError("Matrix order must be positive, got $n"))
    end

    if k < 0
        throw(ArgumentError("Matrix bandwidth must be non-negative, got $k"))
    end

    if num_ccs < 1
        throw(ArgumentError("Minimum number of components must be positive, got $num_ccs"))
    end

    if !(0 < p <= 1)
        throw(ArgumentError("Band density must be in (0, 1], got $p"))
    end

    if n < k + num_ccs
        throw(ArgumentError("Need `n ≥ k + num_ccs`, got n=$n, k=$k, and num_css=$num_ccs"))
    end

    cc_sizes = fill(1, num_ccs)
    leftover = n - num_ccs

    if leftover > 0
        cc_sizes_prop = rand(num_ccs)
        cc_sizes .+= Int.(round.(cc_sizes_prop ./ sum(cc_sizes_prop) .* leftover))
    end

    cc_sizes[1] += n - sum(cc_sizes) # Ensure the sum is exactly `n`

    if cc_sizes[1] <= k
        diff = k - cc_sizes[1] + 1
        cc_sizes[1] += diff
        i = 2

        while diff > 0
            if cc_sizes[i] > 1
                shrink = min(diff, cc_sizes[i] - 1)
                cc_sizes[i] -= shrink
                diff -= shrink
            end

            i += 1
        end
    end

    components = Vector{SparseMatrixCSC{Float64,Int64}}(undef, num_ccs);
    # Ensure that at least one component (arbitrarily, the first) has bandwidth exactly `k`
    components[1] = sparse(random_banded_matrix(cc_sizes[1], k; p=p));

    for i in 2:num_ccs # Some components may themselves be disconnected
        cc_size = cc_sizes[i]
        cc_band = rand(0:min(k, cc_size - 1));
        components[i] = sparse(random_banded_matrix(cc_size, cc_band; p=p));
    end

    #= This matrix will have bandwidth precisely `k` with multiple connected components,
    acting as a stress test for the (reverse) Cuthill–McKee algorithm. =#
    return blockdiag(components...)
end
