# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestUtils

Test suite for the `random_banded_matrix` function.

[TODO: Elaborate on why this is the only one?]
"""
module TestUtils

using MatrixBandwidth
using Random
using Test

const MAX_ORDER = 20

# Assuming that `A` is square and `-n < k < n`
function kth_diagonal(A::Matrix{Float64}, k::Int)
    return map(i -> A[i - min(k, 0), i + max(k, 0)], 1:(size(A, 1) - abs(k)))
end

kth_diagonal_is_empty(A::Matrix{Float64}, k::Int) = all(iszero, kth_diagonal(A, k))

kth_diagonal_is_full(A::Matrix{Float64}, k::Int) = all(!iszero, kth_diagonal(A, k))

function kth_diagonal_has_one_nonzero(A::Matrix{Float64}, k::Int)
    return count(!iszero, kth_diagonal(A, k)) == 1
end

@testset "Random banded matrices – Default density" begin
    for n in 1:MAX_ORDER, k in 0:(n - 1)
        A = random_banded_matrix(n, k)
        @test bandwidth(A) == k
    end
end

@testset "Random banded matrices – Random densities" begin
    for n in 1:MAX_ORDER, k in 0:(n - 1)
        A = random_banded_matrix(n, k; p=rand())
        @test bandwidth(A) == k
    end
end

@testset "Random banded matrices – With RNGs" begin
    rng = MersenneTwister(228)

    for n in 1:MAX_ORDER, k in 0:(n - 1)
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

@testset "Random banded matrices – Sparse bands" begin
    # A sufficiently small probability to ensure that every `rand()` call falls below it
    epsilon = 1 / typemax(UInt128)

    for n in 1:MAX_ORDER, k in 0:(n - 1)
        A = random_banded_matrix(n, k; p=epsilon)
        @test all(map(d -> kth_diagonal_has_one_nonzero(A, d), (-k):k))
        @test all(map(d -> kth_diagonal_is_empty(A, d), (1 - n):(-k - 1))) # Superdiagonals
        @test all(map(d -> kth_diagonal_is_empty(A, d), (k + 1):(n - 1))) # Subdiagonals
    end
end

@testset "Random banded matrices – Full bands" begin
    for n in 1:MAX_ORDER, k in 0:(n - 1)
        A = random_banded_matrix(n, k; p=1)
        @test all(map(d -> kth_diagonal_is_full(A, d), (-k):k))
        @test all(map(d -> kth_diagonal_is_empty(A, d), (1 - n):(-k - 1))) # Superdiagonals
        @test all(map(d -> kth_diagonal_is_empty(A, d), (k + 1):(n - 1))) # Subdiagonals
    end
end

end
