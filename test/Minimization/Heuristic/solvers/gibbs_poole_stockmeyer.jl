# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

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

const MAX_ORDER = 150
const MAX_BAND = 15
const MAX_DENSITY = 0.5
const MAX_CCS = 3
const TEST_PROB = 0.3

include("../test_utils.jl")

#= Check that the algorithm is generally effective on random banded matrices whose rows and
columns have been shuffled. We vary the matrix size `n`, the number of connected components
`num_ccs` (as a stress test for Gibbs–Poole–Stockmeyer, which processes each component
independently), and the original bandwidth `k`. We also vary the density of the off-diagonal
bands (pre-shuffling). The default Hou, Liu, and Zhu (2024) node finder is used .=#
@testset "GPS solver (default RCM++ node finder) – Random matrices (n ≤ $MAX_ORDER)" begin
    Random.seed!(7874336)

    sum_mult_errors = 0.0
    num_cases = 0

    for n in 2:MAX_ORDER
        if rand() < TEST_PROB # Randomly skip some orders to reduce test time
            num_ccs = rand(1:min(MAX_CCS, n - 1))
            k = rand(1:min(MAX_BAND, n - num_ccs))
            p = MAX_DENSITY * rand()

            A = random_banded_discon_matrix(n, k, num_ccs, p)
            perm = randperm(n)
            A = A[perm, perm]
            res = minimize_bandwidth(A, GibbsPooleStockmeyer())

            #= Confirm that the ratio of the minimized bandwidth to the original is not too
            large. (In very rare cases, it is possible that the original ordering of `A` is
            not optimal and that Gibbs–Poole–Stockmeyer produces an even better one.) =#
            mult_error = res.bandwidth / k
            @test mult_error < 2.5

            sum_mult_errors += mult_error
            num_cases += 1

            # Confirm that the ordering computed indeed yields the alleged bandwidth
            @test bandwidth(A[res.ordering, res.ordering]) == res.bandwidth
        end
    end

    #= We keep track of the average ratio of the minimized bandwidth to the original
    bandwidth, which should be less than 1.5 if the algorithm is working properly. =#
    @test sum_mult_errors / num_cases < 1.25
end

#= Here we specify the node finder, using the traditional George and Liu (1979) algorithm
(upon which our default node finder is an improvement). =#
@testset "GPS solver (George–Liu node finder) – Random matrices (n ≤ $MAX_ORDER)" begin
    Random.seed!(101304)
    gl_node_finder = Minimization.Heuristic.pseudo_peripheral_node_finder

    sum_mult_errors = 0.0
    num_cases = 0

    for n in 2:MAX_ORDER
        if rand() < TEST_PROB # Randomly skip some orders to reduce test time
            num_ccs = rand(1:min(MAX_CCS, n - 1))
            k = rand(1:min(MAX_BAND, n - num_ccs))
            p = MAX_DENSITY * rand()

            A = random_banded_discon_matrix(n, k, num_ccs, p)
            perm = randperm(n)
            A = A[perm, perm]
            res = minimize_bandwidth(A, GibbsPooleStockmeyer(gl_node_finder))

            mult_error = res.bandwidth / k
            @test mult_error < 2.5

            sum_mult_errors += mult_error
            num_cases += 1

            @test bandwidth(A[res.ordering, res.ordering]) == res.bandwidth
        end
    end

    @test sum_mult_errors / num_cases < 1.25
end

#= Now a (far inferior) naive node finder is used, simply selecting the node with the lowest
vertex degree per connected component as the starting point. We see that
Gibbs–Poole–Stockmeyer continues to prove effective (although less so) nonetheless. =#
@testset "GPS solver (naive node finder) – Random matrices (n ≤ $MAX_ORDER)" begin
    Random.seed!(48664)
    naive_node_finder =
        A::AbstractMatrix{Bool} -> argmin(i -> sum(view(A, :, i)), axes(A, 1))

    sum_mult_errors = 0.0
    num_cases = 0

    for n in 2:MAX_ORDER
        if rand() < TEST_PROB # Randomly skip some orders to reduce test time
            num_ccs = rand(1:min(MAX_CCS, n - 1))
            k = rand(1:min(MAX_BAND, n - num_ccs))
            p = MAX_DENSITY * rand()

            A = random_banded_discon_matrix(n, k, num_ccs, p)
            perm = randperm(n)
            A = A[perm, perm]
            res = minimize_bandwidth(A, GibbsPooleStockmeyer(naive_node_finder))

            mult_error = res.bandwidth / k
            @test mult_error < 2.5

            sum_mult_errors += mult_error
            num_cases += 1

            @test bandwidth(A[res.ordering, res.ordering]) == res.bandwidth
        end
    end

    @test sum_mult_errors / num_cases < 1.25
end

end
