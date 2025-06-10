# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    _assert_valid_node_selector(selector) -> Nothing

TODO: Write here
"""
function _assert_valid_node_selector(selector::Function)
    input = [false true; true false]

    if !hasmethod(selector, Tuple{AbstractMatrix{Bool}})
        throw(ArgumentError("`selector` must take an an AbstractMatrix{Bool} as input"))
    end

    try
        output = selector(input)
        if !(output isa Integer)
            throw(ArgumentError("`selector` must return an Integer"))
        end
    catch
        throw(
            ArgumentError(
                "`selector` throws an error when called on our test input (K₂'s adjacency)"
            ),
        )
    end

    return nothing
end

"""
    pseudo_peripheral_node(A::AbstractMatrix{Bool}) -> Int

TODO: Write here
"""
function pseudo_peripheral_node(A::AbstractMatrix{Bool})
    n = size(A, 1)

    if n == 1
        return 1
    end

    function _farthest!(distances::Vector{Int}, queue::Queue{Int}, v::Int)
        fill!(distances, typemax(Int))
        distances[v] = 0
        empty!(queue)
        enqueue!(queue, v)

        while !isempty(queue)
            curr = dequeue!(queue)
            neighbors = findall(A[curr, :])

            for neighbor in neighbors
                if distances[neighbor] == typemax(Int)
                    distances[neighbor] = distances[curr] + 1
                    enqueue!(queue, neighbor)
                end
            end
        end

        level = maximum(distances)
        farthest = findall(distances .== level)

        return level, farthest
    end

    function _min_degree_node!(degrees::Vector{Int}, farthest::Vector{Int})
        degrees .= vec(sum(A; dims=1))
        index = argmin(degrees[farthest])
        return farthest[index]
    end

    v = findfirst(!iszero, eachcol(A))
    distances = Vector{Int}(undef, n)
    queue = Queue{Int}()
    level, farthest = _farthest!(distances, queue, v)

    degrees = Vector{Int}(undef, n)
    v_cand = _min_degree_node!(degrees, farthest)
    level_prev = 0

    while level > level_prev
        v, level_prev = v_cand, level
        level, farthest = _farthest!(distances, queue, v)
        v_cand = _min_degree_node!(degrees, farthest)
    end

    return v
end
