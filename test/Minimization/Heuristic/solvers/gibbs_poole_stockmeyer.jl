# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

#= TODO: I feel like this test suite is a bit messy. Let's rewrite it at some point. We
should also add tests for custom node selectors… =#

"""
    TestGibbsPooleStockmeyer

Test suite for the Gibbs–Poole–Stockmeyer matrix bandwidth minimization algorithm.

See `../test_utils.jl` for the definition of `random_banded_discon_matrix`.
"""
module TestGibbsPooleStockmeyer

using MatrixBandwidth
using MatrixBandwidth.Minimization
using Random
using SparseArrays
using Test

const MAX_ORDER = 500
const MAX_BAND = 25
const MAX_DENSITY = 0.5
const MAX_CCS = 10

include("../test_utils.jl")

#= Check that the algorithm is generally effective on random banded matrices whose rows and
columns have been shuffled. We vary the matrix size `n`, the number of connected components
`num_ccs` (as a stress test for Gibbs–Poole–Stockmeyer, which processes each connected
component independently), and the original bandwidth `k`. We also vary the density of the
off-diagonal bands (pre-shuffling). =#
@testset "GPS solver – Random matrices" begin
    Random.seed!(7874336)

    sum_mult_errors = 0.0
    error_count = 0

    for n in 1:MAX_ORDER
        num_ccs = rand(1:min(MAX_CCS, n))
        k = rand(0:min(MAX_BAND, n - num_ccs))
        p = MAX_DENSITY * rand()

        A = random_banded_discon_matrix(n, k, num_ccs, p)
        perm = randperm(n)
        A = A[perm, perm]
        res = minimize_bandwidth(A, GibbsPooleStockmeyer())

        if k > 1 # Avoid division by zero
            #= Confirm that the ratio of the minimized bandwidth to the original is not too
            large. (In very rare cases, it is possible that the original ordering of `A` is
            not optimal and that Gibbs–Poole–Stockmeyer produces an even better one.) =#
            mult_error = res.bandwidth / k
            @test mult_error < 3
            sum_mult_errors += mult_error
            error_count += 1
        end

        # Confirm that the ordering computed indeed yields the alleged bandwidth
        @test bandwidth(A[res.ordering, res.ordering]) == res.bandwidth
    end

    #= We keep track of the average ratio of the minimized bandwidth to the original
    bandwidth, which should be less than 1.5 if the algorithm is working properly. =#
    @test sum_mult_errors / error_count < 1.5
end

end
