# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    HeuristicSolver <: AbstractSolver

Abstract type for all heuristic matrix bandwidth minimization solvers.

Heuristic methods are those which aim to produce near-optimal solutions in a more performant
manner than exact methods. While precise bandwidth minimization is NP-complete, many
heuristic algorithms (such as reverse Cuthill–McKee) run in polynomial time.

Heuristic algorithms differ from metaheuristic ones in that they do not employ higher-level
iterative search frameworks (e.g., stochastic techniques) to survey the global search space
and escape local minima; instead, they rely on straightforward deterministic procedures to
find good solutions in a single pass.
"""
abstract type HeuristicSolver <: AbstractSolver end

_approach(::HeuristicSolver) = :heuristic
