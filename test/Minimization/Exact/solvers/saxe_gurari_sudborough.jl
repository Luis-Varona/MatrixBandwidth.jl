# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestSaxeGurariSudborough

Test suite for the Saxe–Gurari–Sudborough matrix bandwidth minimization algorithm.
"""
module TestSaxeGurariSudborough

using MatrixBandwidth
using MatrixBandwidth.Minimization
using SparseArrays
using Test

const MAX_ORDER = 6
const NUM_ITER = 5

@testset "SGS solver – Brute force verification (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry

        res_bf = minimize_bandwidth(A, BruteForceSearch())
        res_sgs = minimize_bandwidth(A, SaxeGurariSudborough())
        ordering_sgs = res_sgs.ordering

        @test res_bf.bandwidth ==
            res_sgs.bandwidth ==
            bandwidth(A[ordering_sgs, ordering_sgs])
    end
end

end
