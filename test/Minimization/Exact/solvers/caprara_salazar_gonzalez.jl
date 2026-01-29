# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestCapraraSalazarGonzalez

Test suite for the Caprara–Salazar–Gonzalez matrix bandwidth minimization algorithm.
"""
module TestCapraraSalazarGonzalez

using MatrixBandwidth
using MatrixBandwidth.Minimization
using SparseArrays
using Test

const MAX_ORDER = 8
const NUM_ITER = 10

@testset "CSG solver – Brute force verification (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry

        res_bf = minimize_bandwidth(A, BruteForceSearch())
        res_csg = minimize_bandwidth(A, CapraraSalazarGonzalez())
        ordering_csg = res_csg.ordering

        @test res_bf.bandwidth ==
            res_csg.bandwidth ==
            bandwidth(A[ordering_csg, ordering_csg])
    end
end

end
