# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

const DEFAULT_SOLVER = ReverseCuthillMcKee()

"""
    minimize_bandwidth(A::AbstractMatrix) -> BandwidthResult
    minimize_bandwidth(A::AbstractMatrix, solver::AbstractSolver) -> BandwidthResult

TODO: Write here
"""
minimize_bandwidth(A::AbstractMatrix{<:Number}) = minimize_bandwidth(A, DEFAULT_SOLVER)

function minimize_bandwidth(A::AbstractMatrix{<:Number}, solver::AbstractSolver)
    A_bool::AbstractMatrix{Bool} = (!iszero).(A)
    return minimize_bandwidth(A_bool, solver)
end
