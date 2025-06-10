# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestCuthillMcKee

Test suite for the Cuthill–McKee matrix bandwidth minimization algorithm.

See `test/heuristic/utils.jl` for the definitions of `RCM_INPUT` and `RCM_ANSWERS`.
"""
module TestCuthillMcKee

using MatrixBandwidth
using Test

include("utils.jl")

# Check that the algorithm works on a known test case with a symmetric matrix
@testset "CM – Symmetric BitMatrix" begin
    res = minimize_bandwidth(RCM_INPUT::BitMatrix, CuthillMcKee())
    ordering = res.ordering
    #= The answers are given in reverse order, so we cannot trust bandwidth to remain the
    same. Instead, we validate our ordering output only. =#
    @test reverse!(ordering) in last.(RCM_ANSWERS)
end

# Confirm that the solver also accepts asymmetric matrices
@testset "CM – Asymmetric Matrix{Float64}" begin
    A::Matrix{Float64} = rand(10, 10)
    res = minimize_bandwidth(A, CuthillMcKee())
    @test res isa BandwidthResult
end

end
