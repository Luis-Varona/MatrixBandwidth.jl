# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    bi_criteria_node_finder(A) -> Int

Select a pseudo-peripheral node from the graph of `A` based on both eccentricity and width.

This function acts as a node finder for heuristic matrix bandwidth minimization algorithms
such as reverse Cuthill–McKee and Gibbs–Poole–Stockmeyer when applies to connected graphs
(or their adjacency matrices). It heuristically identifies the node "farthest" from the
others in the graph (i.e., a *pseudo-peripheral node*) as a good starting point for the
search process. The heuristic is based on both *eccentricity* (the maximum distance from a
candidate node to any other node in the graph) and *width* (the maximum number of nodes in
any level of the breadth-first search level structure rooted at a candidate node).

It is assumed that `A` is the adjacency matrix of some undirected connected graph;
otherwise, undefined behavior may arise.

# Arguments
- `A::AbstractMatrix{Bool}`: a symmetric boolean matrix with an all-false diagonal,
    representing the adjacency matrix of some undirected connected graph.

# Returns
- `Int`: the index of the pseudo-peripheral node selected from the graph.

# Notes
This algorithm was initially presented in [HLZ24], which in turn is a refinement of the
original, more widely known procedure described in [GL79]. Whereas the [GL79] version used
only eccentricity as a heuristic, the [HLZ24] version uses width as a tiebreaker when
multiple nodes have the same eccentricity, which tends to yield better results.

# References
- [GL79](@cite): A. George and J. W. Liu. *An Implementation of a Pseudoperipheral Node
    Finder*. ACM Transactions on Mathematical Software **5**, 284–95 (1979).
    https://doi.org/10.1145/355841.355845.
- [HLZ24](@cite): J. Hou, H. Liu and S. Zhu. *RCM++:Reverse Cuthill-McKee ordering with
    Bi-Criteria Node Finder* (2024), arXiv:2409.04171 [cs.DS].
    https://arxiv.org/abs/2409.04171
"""
function bi_criteria_node_finder(A::AbstractMatrix{Bool})
    n = size(A, 1)

    if n == 0
        throw(ArgumentError("Matrix order must be positive, got $n"))
    end

    v = findfirst(!iszero, eachcol(A))

    if isnothing(v)
        return 1
    end

    distances = Vector{Int}(undef, n)
    queue = Queue{Int}()
    degrees = Vector{Int}(undef, n)
    record_width = n + 1
    recorded_node = 0

    level, farthest = _farthest!(A, distances, queue, v)
    level_sizes = zeros(Int, n)
    foreach(d -> level_sizes[d + 1] += 1, distances)

    width = maximum(level_sizes)
    v_cand = _min_degree_node!(A, degrees, farthest)
    level_prev = 0

    while level > level_prev
        v = v_cand
        level_prev = level
        level, farthest = _farthest!(A, distances, queue, v)
        fill!(level_sizes, 0)
        foreach(d -> level_sizes[d + 1] += 1, distances)

        width = maximum(level_sizes)
        v_cand = _min_degree_node!(A, degrees, farthest)

        if width < record_width
            record_width = width
            recorded_node = v
        end
    end

    return recorded_node
end

"""
    pseudo_peripheral_node_finder(A) -> Int

Select a pseudo-peripheral node from the graph of `A` based on eccentricity.

This function acts as a node finder for heuristic matrix bandwidth minimization algorithms
such as reverse Cuthill–McKee and Gibbs–Poole–Stockmeyer when applies to connected graphs
(or their adjacency matrices). It heuristically identifies the node "farthest" from the
others in the graph (i.e., a *pseudo-peripheral node*) as a good starting point for the
search process. The heuristic is based solely on *eccentricity* (the maximum distance from a
candidate node to any other node in the graph).

It is assumed that `A` is the adjacency matrix of some undirected connected graph;
otherwise, undefined behavior may arise.

# Arguments
- `A::AbstractMatrix{Bool}`: a symmetric boolean matrix with an all-false diagonal,
    representing the adjacency matrix of some undirected connected graph.

# Returns
- `Int`: the index of the pseudo-peripheral node selected from the graph.

# Notes
This function takes heavy inspiration from the implementation in [Net25], which in turn is
based on the algorithm described in [GL79]. Whereas the [Net25] implementation accepts a
graph object as input and leverages several pre-existing functions in the networkx library,
we repurpose the logic to work directly on adjacency matrices, avoiding reallocation
overhead and an unnecessary dependency on the *Graphs.jl* package.

# References
- [GL79](@cite): A. George and J. W. Liu. *An Implementation of a Pseudoperipheral Node
    Finder*. ACM Transactions on Mathematical Software **5**, 284–95 (1979).
    https://doi.org/10.1145/355841.355845.
- [Net25](@cite): NetworkX Developers. *Source code for networkx.utils.rcm*. NetworkX v3.5
    documentation (2025). Accessed: 2025-06-11.
    https://networkx.org/documentation/stable/_modules/networkx/utils/rcm.html.
"""
function pseudo_peripheral_node_finder(A::AbstractMatrix{Bool})
    n = size(A, 1)

    if n == 0
        throw(ArgumentError("Matrix order must be positive, got $n"))
    end

    v = findfirst(!iszero, eachcol(A))

    if isnothing(v)
        return 1
    end

    distances = Vector{Int}(undef, n)
    queue = Queue{Int}()
    level, farthest = _farthest!(A, distances, queue, v)

    degrees = Vector{Int}(undef, n)
    v_cand = _min_degree_node!(A, degrees, farthest)
    level_prev = 0

    while level > level_prev
        v = v_cand
        level_prev = level
        level, farthest = _farthest!(A, distances, queue, v)
        v_cand = _min_degree_node!(A, degrees, farthest)
    end

    return v
end

function _farthest!(
    A::AbstractMatrix{Bool}, distances::Vector{Int}, queue::Queue{Int}, v::Int
)
    fill!(distances, typemax(Int))
    distances[v] = 0
    empty!(queue)
    enqueue!(queue, v)

    while !isempty(queue)
        curr = dequeue!(queue)
        neighbors = findall(view(A, :, curr))

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

function _min_degree_node!(
    A::AbstractMatrix{Bool}, degrees::Vector{Int}, farthest::Vector{Int}
)
    degrees .= vec(sum(A; dims=1))
    idx = argmin(degrees[farthest])
    return farthest[idx]
end

# Validate that the `node_finder` function is a valid node finder
function _assert_valid_node_finder(node_finder::Function)
    input = [false true; true false]

    if !hasmethod(node_finder, Tuple{AbstractMatrix{Bool}})
        throw(ArgumentError("`node_finder` must take an an AbstractMatrix{Bool} as input"))
    end

    try
        output = node_finder(input)
        if !(output isa Integer)
            throw(ArgumentError("`node_finder` must return an Integer"))
        end
    catch
        throw(
            ArgumentError(
                "`node_finder` throws an error when called on our test input (K₂'s adjacency)",
            ),
        )
    end

    return nothing
end
