# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestCore

Test suite for the core API of the `MatrixBandwidth.jl` package.
"""
module TestCore

using MatrixBandwidth
using MatrixBandwidth.Minimization
using SparseArrays
using Graphs
using Test

const MAX_ORDER1 = 100
const NUM_ITER1 = 100

const MAX_ORDER2 = 8
const NUM_ITER2 = 10

@testset "Testing `bandwidth` (n ≤ $MAX_ORDER1)" begin
    for n in 1:MAX_ORDER1, _ in 1:NUM_ITER1
        density = rand()
        A = sprand(n, n, density)

        k = bandwidth(A)

        @test all(idx -> abs(idx[1] - idx[2]) <= k || A[idx] == 0, CartesianIndices(A))
    end
end

@testset "Testing `bandwidth_lower_bound` (n ≤ $MAX_ORDER2)" begin
    for n in 1:MAX_ORDER2, _ in 1:NUM_ITER2
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Make `A` structurally symmetric

        k = bandwidth_lower_bound(A)
        res = minimize_bandwidth(A, BruteForceSearch())

        @test k <= res.bandwidth
    end
end

@testset "Testing `_floyd_warshall_shortest_paths` (n ≤ $MAX_ORDER1)" begin
    for n in 1:MAX_ORDER1, _ in 1:NUM_ITER1
        density = rand()
        A = sprand(Bool, n, n, density)
        A = A .|| A' # Make `A` structurally symmetric
        g = Graph(A)

        res_matband = MatrixBandwidth._floyd_warshall_shortest_paths(A)
        res_graphs = floyd_warshall_shortest_paths(g).dists
        res_graphs = replace(res_graphs, typemax(Int) => Inf)

        @test res_matband == res_graphs
    end
end

end
