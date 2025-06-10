# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MBPS <: ExactSolver <: AbstractSolver

TODO: Write here
"""
struct MBPS <: ExactSolver
    depth::Int
end

Base.summary(::MBPS) = "Matrix bandwidth by perimeter search"

function _sym_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::MBPS)
    # TODO: Implement
end
