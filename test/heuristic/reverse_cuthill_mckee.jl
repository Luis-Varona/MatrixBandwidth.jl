# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestReverseCuthillMcKee

Test suite for the reverse Cuthill–McKee matrix bandwidth minimization algorithm.

See `test/heuristic/utils.jl` for the definitions of `RCM_INPUT` and `RCM_ANSWERS`.
"""
module TestReverseCuthillMcKee

using MatrixBandwidth
using Random
using SparseArrays
using Test

const RCM_SOLVER = ReverseCuthillMcKee() # Avoid reallocating an identical solver each time
const MAX_ORDER = 200
const MAX_BAND = 10
const MAX_DENSITY = 0.5
const MAX_CCS = 8

include("utils.jl")

#= Check that the algorithm works on a known test case taken from the Boost v1.3.7
documentation. (More details can be found in `test/heuristic/utils.jl`.) =#
@testset "RCM – Known test case" begin
    res = minimize_bandwidth(RCM_INPUT::BitMatrix, RCM_SOLVER)
    answer = (res.bandwidth, res.ordering)
    @test answer in RCM_ANSWERS
end

#= Check that the algorithm is generally effective on random banded matrices whose rows and
columns have been shuffled. We vary the matrix size `n`, the number of connected components
`num_ccs` (as a stress test for reverse Cuthill–McKee, which processes each connected
component independently), and the original bandwidth `k`. We also vary the density of the
off-diagonal bands (pre-shuffling). =#
@testset "RCM – Random matrices" begin
    Random.seed!(7874336)

    sum_mult_errors = 0.0
    error_count = 0

    for n in MAX_CCS:MAX_ORDER, num_ccs in 1:MAX_CCS, k in 0:min(MAX_BAND, n - num_ccs)
        A = random_banded_discon_matrix(n, k, num_ccs, MAX_DENSITY * rand())
        perm = randperm(n)
        A = A[perm, perm]
        res = minimize_bandwidth(A, RCM_SOLVER)

        if k > 1 # Avoid division by zero
            #= Confirm that the ratio of the minimized bandwidth to the original is not too
            large. (In very rare cases, it is possible that the original ordering of `A` is
            not optimal and that reverse Cuthill–McKee produces an even better one.) =#
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
