# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MBID <: ExactSolver <: AbstractSolver

TODO: Write here
"""
struct MBID <: ExactSolver end

Base.summary(::MBID) = "Matrix bandwidth by iterative deepening"

# TODO: Define `minimize_bandwidth` method for `MBID`
