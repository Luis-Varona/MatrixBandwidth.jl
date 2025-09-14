# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestUtils

Test suite for the root-level utility functions of the *MatrixBandwidth.jl* package.
"""
module TestUtils

using MatrixBandwidth
using Graphs
using Random
using SparseArrays
using Test

const MAX_ORDER = 100
const NUM_ITER = 10
const CC_MAX_DENSITY = 0.5
const CC_MAX_NUM_CCS = 5
const N = 100
const P = 0.1

@testset "`random_banded_matrix` – Default density (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, k in 0:(n - 1)
        A = MatrixBandwidth.random_banded_matrix(n, k)
        @test bandwidth(A) == k
        @test MatrixBandwidth.is_structurally_symmetric(A)
    end
end

@testset "`random_banded_matrix` – Random densities (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, k in 0:(n - 1)
        A = MatrixBandwidth.random_banded_matrix(n, k; p=rand())
        @test bandwidth(A) == k
        @test MatrixBandwidth.is_structurally_symmetric(A)
    end
end

@testset "`random_banded_matrix` – With RNGs (n ≤ $MAX_ORDER)" begin
    rng = MersenneTwister(228)

    for n in 1:MAX_ORDER, k in 0:(n - 1)
        A = MatrixBandwidth.random_banded_matrix(n, k; rng=copy(rng))
        B = MatrixBandwidth.random_banded_matrix(n, k; rng=copy(rng))

        p = rand()
        C = MatrixBandwidth.random_banded_matrix(n, k; p=p, rng=copy(rng))
        D = MatrixBandwidth.random_banded_matrix(n, k; p=p, rng=copy(rng))

        @test A == B # Test determinism without a density parameter
        @test C == D # Test determinism with a density parameter

        # `A == B` and `C == D`, so we only need to test `A` and `C`
        @test bandwidth(A) == bandwidth(C) == k
        @test MatrixBandwidth.is_structurally_symmetric(A)
        @test MatrixBandwidth.is_structurally_symmetric(C)
    end
end

@testset "`random_banded_matrix` – Sparse bands (n ≤ $MAX_ORDER)" begin
    #= There is essentially zero chance of any nonzero entries beyond the requisite one per
    superdiagonal and subdiagonal up to the `k`ᵗʰ band. =#
    p = 1 / typemax(UInt128)

    for n in 1:MAX_ORDER, k in 0:(n - 1)
        A = MatrixBandwidth.random_banded_matrix(n, k; p=p)
        @test bandwidth(A) == k

        for i in 1:k
            #= There should be exactly one nonzero entry in each superdiagonal and
            subdiagonal up to the `k`ᵗʰ band. =#
            @test count(!iszero, Iterators.map(j -> A[j, j + i], 1:(n - i))) == 1
            @test count(!iszero, Iterators.map(j -> A[j + i, j], 1:(n - i))) == 1
        end
    end
end

@testset "`random_banded_matrix` – Full bands (n ≤ $MAX_ORDER)" begin
    p = 1

    for n in 1:MAX_ORDER, k in 0:(n - 1)
        A = MatrixBandwidth.random_banded_matrix(n, k; p=p)
        @test bandwidth(A) == k
        # All entries within `k` indices of the diagonal should be nonzero
        @test all(
            idx -> A[idx] != 0,
            Iterators.filter(idx -> abs(idx[1] - idx[2]) <= k, CartesianIndices(A)),
        )
    end
end

@testset "`connected_components` – Singleton" begin
    singleton = Graph(0)
    adj_singleton = Bool.(adjacency_matrix(singleton))
    ccs_singleton = MatrixBandwidth.connected_components(adj_singleton)

    @test isempty(ccs_singleton)
end

@testset "`connected_components` – Empty graphs (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER
        e_n = Graph(n)
        adj_e_n = Bool.(adjacency_matrix(e_n))
        ccs_e_n = MatrixBandwidth.connected_components(adj_e_n)

        @test length(ccs_e_n) == n
        @test all(cc -> length(cc) == 1, ccs_e_n)
        @test sort!(vcat(ccs_e_n...)) == 1:n
    end
end

@testset "`connected_components` – Complete graphs (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER
        k_n = complete_graph(n)
        adj_k_n = Bool.(adjacency_matrix(k_n))
        ccs_k_n = MatrixBandwidth.connected_components(adj_k_n)

        @test length(ccs_k_n) == 1
        @test sort!(vcat(ccs_k_n...)) == 1:n
    end
end

@testset "`connected_components` – Random graphs (n ≤ $MAX_ORDER)" begin
    Random.seed!(122105)

    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        g = erdos_renyi(n, CC_MAX_DENSITY * rand())
        adj_g = Bool.(adjacency_matrix(g))
        ccs_g = MatrixBandwidth.connected_components(adj_g)
        trusted_ccs_g = connected_components(g) # From `Graphs.jl`

        @test sort!(map(sort!, ccs_g)) == sort!(map(sort!, trusted_ccs_g))
    end
end

max_n_next_test = MAX_ORDER * CC_MAX_NUM_CCS

@testset "`connected_components` – Random disconnected graphs (n ≤ $max_n_next_test)" begin
    Random.seed!(040307)

    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        num_ccs = rand(1:CC_MAX_NUM_CCS)
        gs = map(_ -> erdos_renyi(n, CC_MAX_DENSITY * rand()), 1:num_ccs)
        g = reduce(blockdiag, gs)
        ccs_g = MatrixBandwidth.connected_components(Bool.(adjacency_matrix(g)))
        trusted_ccs_g = connected_components(g) # From `Graphs.jl`

        @test sort!(map(sort!, ccs_g)) == sort!(map(sort!, trusted_ccs_g))
    end
end

@testset "`floyd_warshall_shortest_paths` (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        density = rand()
        A = sprand(Bool, n, n, density)
        A = A .|| A' # Ensure structural symmetry
        g = Graph(A)

        res_matband = MatrixBandwidth.floyd_warshall_shortest_paths(A)
        res_graphs = Graphs.floyd_warshall_shortest_paths(g).dists
        res_graphs = replace(res_graphs, typemax(Int) => Inf)

        @test res_matband == res_graphs
    end
end

@testset "`is_structurally_symmetric`" begin
    rng = MersenneTwister(787)
    A = sprand(rng, N, N, P)
    A = A + A' # Ensure structural symmetry

    B = copy(A)
    offdiag_nonzero_idx = first(
        Iterators.filter(idx -> idx[1] != idx[2] && B[idx] != 0, CartesianIndices(B))
    )
    B[offdiag_nonzero_idx] = 0 # Break structural symmetry

    @test MatrixBandwidth.is_structurally_symmetric(A)
    @test !MatrixBandwidth.is_structurally_symmetric(B)
end

@testset "`offdiag_nz_support`" begin
    rng = MersenneTwister(887853)
    A = sprand(rng, N, N, P)

    offdiag_nonzero_idxs = MatrixBandwidth.offdiag_nz_support(A)
    offdiag_zero_idxs = (!).(offdiag_nonzero_idxs)
    foreach(i -> offdiag_zero_idxs[i, i] = false, 1:N)

    @test all(!iszero, A[offdiag_nonzero_idxs])
    @test iszero(A[offdiag_zero_idxs])
end

@testset "`find_direct_subtype`" begin
    abstract type Parent end

    abstract type Child1 <: Parent end
    abstract type Grandchild1 <: Child1 end
    struct Grandchild2 <: Child1 end

    abstract type Child2 <: Parent end
    struct Child3 <: Parent end

    @test MatrixBandwidth.find_direct_subtype(Parent, Child1) === Child1
    @test MatrixBandwidth.find_direct_subtype(Parent, Grandchild1) === Child1
    @test MatrixBandwidth.find_direct_subtype(Parent, Grandchild2) === Child1
    @test MatrixBandwidth.find_direct_subtype(Parent, Child2) === Child2
    @test MatrixBandwidth.find_direct_subtype(Parent, Child3) === Child3
end

end
