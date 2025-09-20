# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestCapraraSalazarGonzalez

Test suite for the Caprara–Salazar–Gonzalez matrix bandwidth recognition algorithm.
"""
module TestCapraraSalazarGonzalez

using MatrixBandwidth
using MatrixBandwidth.Recognition
using SparseArrays
using Test

const MAX_ORDER = 7
const NUM_ITER = 10

@testset "CSG decider – Bandwidth ≤ k (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry

        b = minimize_bandwidth(A, Minimization.BruteForceSearch()).bandwidth
        k = rand(b:(n - 1))

        res = has_bandwidth_k_ordering(A, k, CapraraSalazarGonzalez())
        ordering = res.ordering

        @test res.has_ordering
        @test bandwidth(A[ordering, ordering]) <= k
    end
end

@testset "CSG decider – Bandwidth > k (n ≤ $MAX_ORDER)" begin
    for n in 2:MAX_ORDER, i in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry

        while bandwidth(A) == 0
            density = rand()
            A = sprand(n, n, density)
            A = A + A' # Ensure structural symmetry
        end

        b = minimize_bandwidth(A, Minimization.BruteForceSearch()).bandwidth
        k = rand(0:(b - 1))

        res = has_bandwidth_k_ordering(A, k, CapraraSalazarGonzalez())

        @test !res.has_ordering
        @test isnothing(res.ordering)
    end
end

end
