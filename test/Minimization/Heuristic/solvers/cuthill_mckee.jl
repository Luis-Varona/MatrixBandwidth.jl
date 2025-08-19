# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestCuthillMcKee

Test suite for the Cuthill–McKee and reverse Cuthill–McKee matrix bandwidth minimization
algorithms.

See `../test_utils.jl` for the definitions of `RCM_INPUT`, `RCM_ANSWERS`, and
`random_banded_discon_matrix`.
"""
module TestCuthillMcKee

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

#= Check that the algorithm works on a known test case taken from the Boost v1.3.7
documentation. (More details can be found in `test/heuristic/utils.jl`.) =#
@testset "CM solver – Known test case" begin
    res = minimize_bandwidth(RCM_INPUT::BitMatrix, CuthillMcKee())
    ordering = res.ordering
    #= The answers are given in reverse order, so we cannot trust bandwidth to remain the
    same. Instead, we validate our ordering output only. =#
    @test reverse!(ordering) in last.(RCM_ANSWERS)
end

@testset "RCM solver – Known test case" begin
    res = minimize_bandwidth(RCM_INPUT::BitMatrix, ReverseCuthillMcKee())
    answer = (res.bandwidth, res.ordering)
    @test answer in RCM_ANSWERS
end

#= NOTE: We only bother testing RCM henceforth, since it is already a wrapper around CM that
induces the same bandwidth in the end. =#

#= Check that the algorithm is generally effective on random banded matrices whose rows and
columns have been shuffled. We vary the matrix size `n`, the number of connected components
`num_ccs` (as a stress test for reverse Cuthill–McKee, which processes each component
independently), and the original bandwidth `k`. We also vary the density of the off-diagonal
bands (pre-shuffling). The default Hou, Liu, and Zhu (2024) node finder is used .=#
@testset "RCM solver (default RCM++ node finder) – Random matrices (n ≤ $MAX_ORDER)" begin
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
            res = minimize_bandwidth(A, ReverseCuthillMcKee())

            #= Confirm that the ratio of the minimized bandwidth to the original is not too
            large. (In very rare cases, it is possible that the original ordering of `A` is not
            optimal and that reverse Cuthill–McKee produces an even better one.) =#
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
@testset "RCM solver (George–Liu node finder) – Random matrices (n ≤ $MAX_ORDER)" begin
    Random.seed!(56674)
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
            res = minimize_bandwidth(A, ReverseCuthillMcKee(gl_node_finder))

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
@testset "RCM solver (naive node finder) – Random matrices (n ≤ $MAX_ORDER)" begin
    Random.seed!(84664)
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
