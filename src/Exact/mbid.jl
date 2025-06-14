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

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::MBID)
    n = size(A, 1)
    adj_lists = map(node -> findall(A[:, node]), 1:n)

    function _add_node!(
        ordering_buf::Vector{Int},
        unselected::Set{Int},
        adj_list::Set{Int},
        bandwidth::Int,
        depth::Int,
    )
        if depth == n
            ordering = ordering_buf
        else
            if depth == 0
                starting_node = minimum(unselected)
            end

            ordering = nothing
            res = iterate(unselected)

            while (isnothing(ordering) && !isnothing(res))
                (node, state) = res

                if depth > 0 || node == starting_node
                    ordering_buf[depth + 1] = node
                    unselected_new = setdiff(unselected, [node])
                    adj_list_new = union(adj_list, adj_lists[node])
                    adj_list_new = intersect(adj_list_new, unselected_new)

                    if _is_compatible(A, ordering_buf, adj_list_new, bandwidth, depth)
                        ordering = _add_node!(
                            ordering_buf, unselected_new, adj_list_new, bandwidth, depth + 1
                        )
                    end
                end

                res = iterate(unselected, state)
            end
        end

        return ordering
    end

    ordering = nothing
    ordering_buf = Vector{Int}(undef, n)
    bandwidth = ceil(Int, maximum(length, adj_lists) / 2)

    while isnothing(ordering)
        unselected = Set(1:n)
        adj_list = Set{Int}()
        ordering = _add_node!(ordering_buf, unselected, adj_list, bandwidth, 0)
        bandwidth += 1
    end

    return ordering
end

function _is_compatible(
    A::AbstractMatrix{Bool},
    ordering::Vector{Int},
    adj_list::Set{Int},
    bandwidth::Int,
    depth::Int,
)
    # TODO: Implement
    return true # Just a placeholder
end
