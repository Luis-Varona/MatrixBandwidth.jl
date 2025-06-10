# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    minimize_bandwidth(A, solver=ReverseCuthillMcKee()) -> BandwidthResult

TODO: Write here
"""
function minimize_bandwidth(
    A::AbstractMatrix{T}, solver::AbstractSolver=DEFAULT_SOLVER
) where {T<:Number}
    if !allequal(size(A))
        throw(ArgumentError("Matrix bandwidth is not defined for non-square matrices"))
    end

    if T !== Bool
        A_bool = (!iszero).(A)
    else
        A_bool = A
    end

    A_copy = copy(A_bool)
    ordering = _sym_minimal_band_ordering(A_copy, solver)
    A_reordered = A_bool[ordering, ordering]
    indices = findall(A_reordered)

    if isempty(indices)
        bandwidth = 0
    else
        bandwidth = maximum(abs(index[1] - index[2]) for index in indices)
    end

    return BandwidthResult(A, bandwidth, ordering, solver)
end

"""
    _sym_minimal_bandwidh_ordering(A::AbstractMatrix{Bool}, solver::AbstractSolver)
        -> Vector{Int}

TODO: Write here
"""
function _sym_minimal_band_ordering(::AbstractMatrix{Bool}, ::T) where {T<:AbstractSolver}
    S = supertype(T)

    if S === AbstractSolver
        subtype = T
    else
        subtype = S
    end

    throw(NotImplementedError(approach, :solver, subtype, AbstractSolver))
end
