# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestMinimization

Test suite for the `Minimization` submodule of the MatrixBandwidth.jl package.
"""
module TestMinimization

const TEST_GROUPS = []
const NESTED_TEST_SUITES = [
    "Exact/Exact.jl", "Heuristic/Heuristic.jl", "Metaheuristic/Metaheuristic.jl"
]

for group in TEST_GROUPS
    @info "Testing `Minimization/$group`"
    include("$group.jl")
    println()
end

for suite in NESTED_TEST_SUITES
    include(suite)
end

end
