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
    _assert_matrix_is_square(A)

    if T !== Bool
        A_bool = (!iszero).(A)
    else
        A_bool = copy(A)
    end

    ordering = _bool_minimal_band_ordering(A_bool, solver)
    A_reordered = A_bool[ordering, ordering]
    bandwidth = bandwidth_unpermuted(A_reordered)

    return BandwidthResult(A, bandwidth, ordering, solver)
end

"""
    bandwidth_unpermuted(A) -> Int

TODO: Write here
"""
function bandwidth_unpermuted(A::AbstractMatrix{T}) where {T<:Number}
    _assert_matrix_is_square(A)

    if T !== Bool
        A_bool = (!iszero).(A)
    else
        A_bool = A
    end

    indices = findall(A_bool)

    if isempty(indices)
        bandwidth = 0
    else
        bandwidth = maximum(abs(index[1] - index[2]) for index in indices)
    end

    return bandwidth
end

#= Compute a minimal bandwidth ordering for a preprocessed `AbstractMatrix{Bool}`.
Restricting entries to booleans can improve performance via cache optimizations, bitwise
operations, etc. Each concrete `AbstractSolver` subtype must implement its own
`_bool_minimal_band_ordering` method to define the corresponding algorithm logic. =#
function _bool_minimal_band_ordering(::AbstractMatrix{Bool}, ::T) where {T<:AbstractSolver}
    S = supertype(T)

    if S === AbstractSolver
        subtype = T
    else
        subtype = S
    end

    throw(NotImplementedError(approach, :solver, subtype, AbstractSolver))
end
