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
