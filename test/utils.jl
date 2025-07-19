# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestUtils

[TODO: Write here.]
"""
module TestUtils

using MatrixBandwidth
using Random
using Test

const MAX_ORDER = 20

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

# TODO: Add tests for `random_banded_matrix` with sparse bands (`p = 1 / typemax(UInt128)`)

# TODO: Add tests for `random_banded_matrix` with full bands (`p = 1`)

# TODO: Add tests for `_find_direct_subtype`

# TODO: Add tests for `_is_structurally_symmetric`

# TODO: Add tests for `_offdiag_nonzero_support`

end
