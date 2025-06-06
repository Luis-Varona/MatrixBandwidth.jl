# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    GRASP <: MetaheuristicSolver <: AbstractSolver

TODO: Write here
"""
struct GRASP <: MetaheuristicSolver end

Base.summary(::GRASP) = "Greedy randomized adaptive search procedure"

# TODO: Define `minimize_bandwidth` method for `GRASP`
