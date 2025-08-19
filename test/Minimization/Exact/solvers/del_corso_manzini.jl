# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestDelCorsoManzini

Test suite for the Del Corso–Manzini matrix bandwidth minimization algorithm, with and
without the incorporation of perimeter search.
"""
module TestDelCorsoManzini

using MatrixBandwidth
using MatrixBandwidth.Minimization
using SparseArrays
using Test

const MAX_ORDER = 6
const NUM_ITER = 50

@testset "DCM solver – Brute force verification (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry

        res_bf = minimize_bandwidth(A, BruteForceSearch())
        res_dcm = minimize_bandwidth(A, DelCorsoManzini())
        ordering_dcm = res_dcm.ordering

        @test res_bf.bandwidth ==
            res_dcm.bandwidth ==
            bandwidth(A[ordering_dcm, ordering_dcm])
    end
end

@testset "DCM-PS solver (default depth) – Brute force verification (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry

        res_bf = minimize_bandwidth(A, BruteForceSearch())
        res_dcm = minimize_bandwidth(A, DelCorsoManziniWithPS())
        ordering_dcm = res_dcm.ordering

        @test res_bf.bandwidth ==
            res_dcm.bandwidth ==
            bandwidth(A[ordering_dcm, ordering_dcm])
    end
end

@testset "DCM-PS solver (custom depth) – Brute force verification (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry
        depth = rand(1:n)

        res_bf = minimize_bandwidth(A, BruteForceSearch())
        res_dcm = minimize_bandwidth(A, DelCorsoManziniWithPS(depth))
        ordering_dcm = res_dcm.ordering

        @test res_bf.bandwidth ==
            res_dcm.bandwidth ==
            bandwidth(A[ordering_dcm, ordering_dcm])
    end
end

end
