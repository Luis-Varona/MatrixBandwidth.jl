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
    A_sym = _symmetrize(A)
    adj_lists = map(node -> findall(A_sym[:, node]), 1:n)

    function _add_node!(
        ordering_buf::Vector{Int},
        unselected::Set{Int},
        adj_list::Set{Int},
        bandwidth::Int,
        depth::Int,
    )
        if depth == n
            return ordering_buf
        end

        ordering = nothing
        res = iterate(unselected)

        while (isnothing(ordering) && !isnothing(res))
            (node, state) = res
            ordering_buf[depth + 1] = node

            if !any(i -> A_sym[node, ordering_buf[i]] && bandwidth + i <= depth, 1:depth)
                unselected_new = setdiff(unselected, [node])

                if ordering_buf[1] < maximum(unselected_new)
                    adj_list_new = union(adj_list, adj_lists[node])
                    adj_list_new = intersect(adj_list_new, unselected_new)

                    if _is_compatible(A_sym, ordering_buf, adj_list_new, bandwidth, depth)
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
    A_sym::AbstractMatrix{Bool},
    ordering::Vector{Int},
    adj_list::Set{Int},
    bandwidth::Int,
    depth::Int,
)
    if length(adj_list) > bandwidth
        return false
    end

    l = length(adj_list)
    latest_positions = Vector{Int}(undef, l)

    for (i, neighbour) in enumerate(adj_list)
        latest_position = typemax(Int)

        for (j, placed_node) in enumerate(view(ordering, 1:depth))
            if A_sym[neighbour, placed_node]
                latest_position = min(latest_position, bandwidth + j)
            end
        end

        latest_positions[i] = latest_position
    end

    sort!(latest_positions)
    constraints = (1:l) .+ depth

    return all(latest_positions .>= constraints)
end
