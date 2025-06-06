# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SimulatedAnnealing <: MetaheuristicSolver <: AbstractSolver

TODO: Write here
"""
struct SimulatedAnnealing <: MetaheuristicSolver end

Base.summary(::SimulatedAnnealing) = "Simulated annealing"

# TODO: Define `minimize_bandwidth` method for `SimulatedAnnealing`
