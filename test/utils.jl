# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestUtils

Test suite for the root-level utility functions of the `MatrixBandwidth.jl` package.
"""
module TestUtils

using MatrixBandwidth
using Random
using SparseArrays
using Test

const RBM_MAX_ORDER = 20
const N = 100
const P = 0.1

@testset "`random_banded_matrix` – Default density" begin
    for n in 1:RBM_MAX_ORDER, k in 0:(n - 1)
        A = random_banded_matrix(n, k)
        @test bandwidth(A) == k
    end
end

@testset "`random_banded_matrix` – Random densities" begin
    for n in 1:RBM_MAX_ORDER, k in 0:(n - 1)
        A = random_banded_matrix(n, k; p=rand())
        @test bandwidth(A) == k
    end
end

@testset "`random_banded_matrix` – With RNGs" begin
    rng = MersenneTwister(228)

    for n in 1:RBM_MAX_ORDER, k in 0:(n - 1)
        A = random_banded_matrix(n, k; rng=copy(rng))
        B = random_banded_matrix(n, k; rng=copy(rng))

        p = rand()
        C = random_banded_matrix(n, k; p=p, rng=copy(rng))
        D = random_banded_matrix(n, k; p=p, rng=copy(rng))

        @test all(bandwidth.([A, B, C, D]) .== k)
        @test A == B # Test determinism without a density parameter
        @test C == D # Test determinism with a density parameter
    end
end

# TODO: Add tests for `random_banded_matrix` with sparse bands (`p = 1 / typemax(UInt128)`)

# TODO: Add tests for `random_banded_matrix` with full bands (`p = 1`)

@testset "`_find_direct_subtype`" begin
    abstract type Parent end

    abstract type Child1 <: Parent end
    abstract type Grandchild1 <: Child1 end
    struct Grandchild2 <: Child1 end

    abstract type Child2 <: Parent end
    struct Child3 <: Parent end

    @test MatrixBandwidth._find_direct_subtype(Parent, Child1) === Child1
    @test MatrixBandwidth._find_direct_subtype(Parent, Grandchild1) === Child1
    @test MatrixBandwidth._find_direct_subtype(Parent, Grandchild2) === Child1
    @test MatrixBandwidth._find_direct_subtype(Parent, Child2) === Child2
    @test MatrixBandwidth._find_direct_subtype(Parent, Child3) === Child3
end

@testset "`_is_structurally_symmetric`" begin
    rng = MersenneTwister(787)
    A = sprand(rng, N, N, P)
    A = A + A' # Ensure structural symmetry

    B = copy(A)
    offdiag_nonzero_idx = first(
        Iterators.filter(idx -> idx[1] != idx[2] && B[idx] != 0, CartesianIndices(B))
    )
    B[offdiag_nonzero_idx] = 0 # Break structural symmetry

    @test MatrixBandwidth._is_structurally_symmetric(A)
    @test !MatrixBandwidth._is_structurally_symmetric(B)
end

@testset "`_offdiag_nonzero_support`" begin
    rng = MersenneTwister(887853)
    A = sprand(rng, N, N, P)

    offdiag_nonzero_idxs = MatrixBandwidth._offdiag_nonzero_support(A)
    offdiag_zero_idxs = (!).(offdiag_nonzero_idxs)
    foreach(i -> offdiag_zero_idxs[i, i] = false, 1:N)

    @test all(!iszero, A[offdiag_nonzero_idxs])
    @test iszero(A[offdiag_zero_idxs])
end

end
