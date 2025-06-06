# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    GeneticAlgorithm <: MetaheuristicSolver <: AbstractSolver

TODO: Write here
"""
struct GeneticAlgorithm <: MetaheuristicSolver end

Base.summary(::GeneticAlgorithm) = "Genetic algorithm"

# TODO: Define `minimize_bandwidth` method for `GeneticAlgorithm`
