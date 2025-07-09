# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    ExactSolver <: AbstractSolver <: AbstractAlgorithm

Abstract type for all exact matrix bandwidth minimization solvers.

Exact methods are those which guarantee an optimal ordering producing the true minimum
bandwidth of a matrix. Since bandwidth minimization is an NP-complete problem, existing
exact algorithms are, at best, exponential in time complexity—much worse than many
polynomial-time heuristic approaches (e.g., reverse Cuthill–McKee). Such methods, therefore,
are not feasible for large matrices, but they remain useful when precise solutions are
required for small-to-medium-sized inputs (say, up to ``100×100``).
"""
abstract type ExactSolver <: AbstractSolver end

_approach(::ExactSolver) = :exact
