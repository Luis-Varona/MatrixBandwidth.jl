# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestDelCorsoManzini

Test suite for the Del Corso–Manzini matrix bandwidth recognition algorithm.
"""
module TestDelCorsoManzini

using MatrixBandwidth
using MatrixBandwidth.Recognition
using SparseArrays
using Test

const MAX_ORDER = 8
const NUM_ITER = 10

@testset "DCM decider – Bandwidth ≤ k (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry

        b = minimize_bandwidth(A, Minimization.BruteForceSearch()).bandwidth
        k = rand(b:(n - 1))

        res = has_bandwidth_k_ordering(A, k, DelCorsoManzini())
        ordering = res.ordering

        @test res.has_ordering
        @test bandwidth(A[ordering, ordering]) <= k
    end
end

@testset "DCM-PS decider (default depth) – Bandwidth ≤ k (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry

        b = minimize_bandwidth(A, Minimization.BruteForceSearch()).bandwidth
        k = rand(b:(n - 1))

        res = has_bandwidth_k_ordering(A, k, DelCorsoManziniWithPS())
        ordering = res.ordering

        @test res.has_ordering
        @test bandwidth(A[ordering, ordering]) <= k
    end
end

@testset "DCM-PS decider (custom depth) – Bandwidth ≤ k (n ≤ $MAX_ORDER)" begin
    for n in 1:MAX_ORDER, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry

        b = minimize_bandwidth(A, Minimization.BruteForceSearch()).bandwidth
        k = rand(b:(n - 1))
        depth = rand(1:n)

        res = has_bandwidth_k_ordering(A, k, DelCorsoManziniWithPS(depth))
        ordering = res.ordering

        @test res.has_ordering
        @test bandwidth(A[ordering, ordering]) <= k
    end
end

@testset "DCM decider – Bandwidth > k (n ≤ $MAX_ORDER)" begin
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

        res = has_bandwidth_k_ordering(A, k, DelCorsoManzini())

        @test !res.has_ordering
        @test isnothing(res.ordering)
    end
end

@testset "DCM-PS decider (default depth) – Bandwidth > k (n ≤ $MAX_ORDER)" begin
    for n in 2:MAX_ORDER, _ in 1:NUM_ITER
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

        res = has_bandwidth_k_ordering(A, k, DelCorsoManziniWithPS())

        @test !res.has_ordering
        @test isnothing(res.ordering)
    end
end

@testset "DCM-PS decider (custom depth) – Bandwidth > k (n ≤ $MAX_ORDER)" begin
    for n in 2:MAX_ORDER, _ in 1:NUM_ITER
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
        depth = rand(1:n)

        res = has_bandwidth_k_ordering(A, k, DelCorsoManziniWithPS(depth))

        @test !res.has_ordering
        @test isnothing(res.ordering)
    end
end

end
