# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    DelCorsoManzini <: AbstractDecider

TODO: Write here
"""
struct DelCorsoManzini <: AbstractDecider end

Base.summary(::DelCorsoManzini) = "Del Corso–Manzini algorithm"

function _bool_bandwidth_k_ordering(A::AbstractMatrix{Bool}, k::Int, ::DelCorsoManzini)
    n = size(A, 1)

    ordering_buf = Vector{Int}(undef, n)
    A_sym = _symmetrize(A)
    adj_lists = map(node -> findall(A_sym[:, node]), 1:n)

    unselected = Set(1:n)
    adj_list = Set{Int}()
    depth = 0

    return _add_node!(ordering_buf, A_sym, k, adj_lists, unselected, adj_list, depth)
end

# TODO: Remove `depth` parameter in favor of computing `n - length(unselected)` each time?
function _add_node!(
    ordering_buf::Vector{Int},
    A_sym::AbstractMatrix{Bool},
    k::Int,
    adj_lists::Vector{Vector{Int}},
    unselected::Set{Int},
    adj_list::Set{Int},
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

        if !any(i -> A_sym[node, ordering_buf[i]] && k + i <= depth, 1:depth)
            unselected_new = setdiff(unselected, [node])

            if ordering_buf[1] < maximum(unselected_new)
                adj_list_new = union(adj_list, adj_lists[node])
                adj_list_new = intersect(adj_list_new, unselected_new)

                if _is_compatible(A_sym, ordering_buf, adj_list_new, k, depth)
                    ordering = _add_node!(
                        ordering_buf,
                        A_sym,
                        k,
                        adj_lists,
                        unselected_new,
                        adj_list_new,
                        depth + 1,
                    )
                end
            end

            res = iterate(unselected, state)
        end
    end

    return ordering
end
