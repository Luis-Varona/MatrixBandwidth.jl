# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MetaheuristicSolver <: AbstractSolver <: AbstractAlgorithm

Abstract type for all metaheuristic matrix bandwidth minimization solvers.

Metaheuristic methods are those which employ higher-level iterative search frameworks such
as stochastic techniques or nature-inspired processes to survey the global search space and
escape local minima. Unlike heuristic methods—which follow fixed deterministic
procedures—metaheuristics adaptively refine candidate solutions over multiple iterations.
Although metaheuristic approaches are often slower than heuristic ones (but certainly still
faster than exact ones), they shine in complex cases where the latter may get trapped in
poor-quality local minima.

# Supertype Hierarchy
`MetaheuristicSolver` <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)
"""
abstract type MetaheuristicSolver <: AbstractSolver end

Minimization._approach(::MetaheuristicSolver) = :metaheuristic
