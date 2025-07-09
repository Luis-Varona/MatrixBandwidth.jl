# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    DelCorsoManzini <: AbstractDecider <: AbstractAlgorithm

TODO: Write here
"""
struct DelCorsoManzini <: AbstractDecider end

Base.summary(::DelCorsoManzini) = "Del Corso–Manzini"

_requires_symmetry(::DelCorsoManzini) = true

function _bool_bandwidth_k_ordering(A::AbstractMatrix{Bool}, k::Int, ::DelCorsoManzini)
    n = size(A, 1)

    ordering_buf = Vector{Int}(undef, n)
    adj_lists = map(node -> findall(A[:, node]), 1:n)

    unselected = Set(1:n)
    adj_list = Set{Int}()
    depth = 0

    return _add_node!(
        ordering_buf, A, k, adj_lists, unselected, adj_list, depth, DelCorsoManzini()
    )
end

function _add_node!(
    ordering_buf::Vector{Int},
    A::AbstractMatrix{Bool},
    k::Int,
    adj_lists::Vector{Vector{Int}},
    unselected::Set{Int},
    adj_list::Set{Int},
    depth::Int,
    decider::DelCorsoManzini,
)
    if depth == size(A, 1)
        return ordering_buf
    end

    ordering = nothing
    res = iterate(unselected)

    while (isnothing(ordering) && !isnothing(res))
        (node, state) = res
        ordering_buf[depth + 1] = node

        if !any(i -> A[node, ordering_buf[i]] && k + i <= depth, 1:depth)
            unselected_new = setdiff(unselected, [node])

            if isempty(unselected_new)
                max_label = node
            else
                max_label = maximum(unselected_new)
            end

            if ordering_buf[1] < max_label
                adj_list_new = union(adj_list, adj_lists[node])
                adj_list_new = intersect(adj_list_new, unselected_new)

                if _is_compatible(A, ordering_buf, adj_list_new, k, depth)
                    ordering = _add_node!(
                        ordering_buf,
                        A,
                        k,
                        adj_lists,
                        unselected_new,
                        adj_list_new,
                        depth + 1,
                        decider,
                    )
                end
            end
        end

        res = iterate(unselected, state)
    end

    return ordering
end

function _is_compatible(
    A::AbstractMatrix{Bool}, ordering::Vector{Int}, adj_list::Set{Int}, k::Int, depth::Int
)
    if length(adj_list) > k
        return false
    end

    l = length(adj_list)
    latest_positions = Vector{Int}(undef, l)

    for (i, neighbor) in enumerate(adj_list)
        latest_position = typemax(Int)

        for (j, placed_node) in enumerate(view(ordering, 1:depth))
            if A[neighbor, placed_node]
                latest_position = min(latest_position, k + j)
            end
        end

        latest_positions[i] = latest_position
    end

    sort!(latest_positions)
    constraints = (1:l) .+ depth

    return all(latest_positions .>= constraints)
end
