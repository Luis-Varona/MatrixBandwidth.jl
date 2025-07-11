# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestMetaheuristic

Test suite for the `Minimization.Metaheuristic` submodule of the `MatrixBandwidth.jl`
package.
"""
module TestMetaheuristic

include("solvers/grasp.jl")
include("solvers/simulated_annealing.jl")
include("solvers/genetic_algorithm.jl")

end
