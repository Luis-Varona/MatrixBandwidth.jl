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
using Test

include("utils.jl")

# Check that the algorithm works on a known test case with a symmetric matrix
@testset "RCM – Symmetric BitMatrix" begin
    res = minimize_bandwidth(RCM_INPUT::BitMatrix, ReverseCuthillMcKee())
    answer = (res.bandwidth, res.ordering)
    @test answer in RCM_ANSWERS
end

# Confirm that the solver also accepts asymmetric matrices
@testset "RCM – Asymmetric Matrix{Float64}" begin
    A::Matrix{Float64} = rand(10, 10)
    res = minimize_bandwidth(A, ReverseCuthillMcKee())
    @test res isa BandwidthResult
end

end
