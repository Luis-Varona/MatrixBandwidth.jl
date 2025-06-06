# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

const DEFAULT_SOLVER = ReverseCuthillMcKee()

"""
    minimize_bandwidth(A, solver=ReverseCuthillMcKee()) -> BandwidthResult

TODO: Write here
"""
function minimize_bandwidth(
    A::AbstractMatrix{<:Bool}, solver::AbstractSolver=DEFAULT_SOLVER
)
    if !allequal(size(A))
        throw(ArgumentError("Matrix bandwidth is not defined for non-square matrices"))
    end

    A_copy = copy(A)
    return _minimize_bandwidth_safe(A_copy, solver)
end

function minimize_bandwidth(
    A::AbstractMatrix{<:Number}, solver::AbstractSolver=DEFAULT_SOLVER
)
    A_bool::AbstractMatrix{Bool} = (!iszero).(A)
    return minimize_bandwidth(A_bool, solver)
end
