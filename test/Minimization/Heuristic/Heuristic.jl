# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestHeuristic

Test suite for the `Minimization.Heuristic` submodule of the `MatrixBandwidth.jl` package.
"""
module TestHeuristic

include("solvers/gibbs_poole_stockmeyer.jl")
include("solvers/cuthill_mckee.jl")

end
