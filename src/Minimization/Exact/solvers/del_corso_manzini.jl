# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    DelCorsoManzini <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

TODO: Write here
"""
struct DelCorsoManzini <: ExactSolver end

Base.summary(::DelCorsoManzini) = "Del Corso–Manzini"

_requires_symmetry(::DelCorsoManzini) = true

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::DelCorsoManzini)
    n = size(A, 1)
    ordering = nothing
    ordering_buf = Vector{Int}(undef, n)
    k = bandwidth_lower_bound(A)
    adj_lists = map(node -> findall(A[:, node]), 1:n)
    decider = Recognition.DelCorsoManzini()

    while isnothing(ordering)
        unselected = Set(1:n)
        adj_list = Set{Int}()
        ordering = Recognition._add_node!(
            ordering_buf, A, k, adj_lists, unselected, adj_list, 0, decider
        )
        k += 1
    end

    return ordering
end
