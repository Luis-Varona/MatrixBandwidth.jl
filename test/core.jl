# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestCore

Test suite for the core API of the *MatrixBandwidth.jl* package.
"""
module TestCore

using MatrixBandwidth
using SparseArrays
using Test

const MAX_ORDER1 = 100
const MAX_ORDER2 = 8
const NUM_ITER = 10

@testset "`bandwidth` (n ≤ $MAX_ORDER1)" begin
    for n in 1:MAX_ORDER1, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)

        k = bandwidth(A)

        @test all(idx -> abs(idx[1] - idx[2]) <= k || A[idx] == 0, CartesianIndices(A))
    end
end

@testset "`profile` (n ≤ $MAX_ORDER1)" begin
    for n in 1:MAX_ORDER1, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)

        prof_default = profile(A)
        prof_col = profile(A, dim=:col)
        prof_row = profile(A, dim=:row)

        @test prof_default == prof_col # Should compute column profile by default
        # Only `:col` and `:row` are valid dimensions
        @test_throws ArgumentError profile(A, dim=:nonsense)

        foreach(i -> A[i, i] = 0, 1:n)
        prof_col_zero_diag = profile(A, dim=:col)
        prof_row_zero_diag = profile(A, dim=:row)

        foreach(i -> A[i, i] = 1, 1:n)
        prof_col_one_diag = profile(A, dim=:col)
        prof_row_one_diag = profile(A, dim=:row)

        # Diagonal entries should not affect profile
        @test prof_col == prof_col_zero_diag == prof_col_one_diag
        @test prof_row == prof_row_zero_diag == prof_row_one_diag

        A[:, rand(1:n)] .= 0
        A[rand(1:n), :] .= 0

        # Should run even with zero rows/columns
        @test try
            profile(A);
            true
        catch
            false
        end
    end
end

@testset "`bandwidth_lower_bound` (n ≤ $MAX_ORDER2)" begin
    for n in 1:MAX_ORDER2, _ in 1:NUM_ITER
        density = rand()
        A = sprand(n, n, density)
        A = A + A' # Ensure structural symmetry

        k = bandwidth_lower_bound(A)
        res = minimize_bandwidth(A, Minimization.BruteForceSearch())

        @test k <= res.bandwidth
    end
end

end
