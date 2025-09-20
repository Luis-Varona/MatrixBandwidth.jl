# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestMetaheuristic

Test suite for the `Minimization.Metaheuristic` submodule of the *MatrixBandwidth.jl*
package.
"""
module TestMetaheuristic

const TEST_GROUPS = [
    "solvers/grasp",
    "solvers/psohc",
    "solvers/simulated_annealing",
    "solvers/genetic_algorithm",
    "solvers/ant_colony",
    "solvers/tabu_search",
]

for group in TEST_GROUPS
    @info "Testing `Minimization/Metaheuristic/$group`"
    include("$group.jl")
    println()
end

end
