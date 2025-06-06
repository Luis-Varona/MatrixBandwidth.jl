# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    ReverseCuthillMcKee <: HeuristicSolver <: AbstractSolver

TODO: Write here
"""
struct ReverseCuthillMcKee <: HeuristicSolver end

Base.summary(::ReverseCuthillMcKee) = "Reverse Cuthill–McKee algorithm"

# TODO: Define `minimize_bandwidth` method for `ReverseCuthillMcKee`
