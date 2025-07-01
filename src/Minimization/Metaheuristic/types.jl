# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MetaheuristicSolver <: AbstractSolver

TODO: Write here
"""
abstract type MetaheuristicSolver <: AbstractSolver end

_approach(::MetaheuristicSolver) = :metaheuristic
