# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CuthillMcKee <: HeuristicSolver <: AbstractSolver

TODO: Write here
"""
struct CuthillMcKee <: HeuristicSolver end

Base.summary(::CuthillMcKee) = "Cuthill–McKee algorithm"

function _minimize_bandwidth_safe(A::AbstractMatrix{<:Bool}, ::CuthillMcKee)
    # TODO: Implement
end
