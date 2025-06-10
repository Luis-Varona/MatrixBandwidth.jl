# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CuthillMcKee <: HeuristicSolver <: AbstractSolver

TODO: Write here
"""
struct CuthillMcKee <: HeuristicSolver
    node_selector::Function

    function CuthillMcKee(node_selector::Function=DEFAULT_SELECTOR)
        _assert_valid_node_selector(node_selector)
        return new(node_selector)
    end
end

Base.summary(::CuthillMcKee) = "Cuthill–McKee algorithm"

function _sym_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::CuthillMcKee)
    if A != A'
        A_sym = (!iszero).(A + A') # TODO: Check performance later
    else
        A_sym = A
    end

    n = size(A_sym, 1)
    A_sym[1:(n + 1):end] .= false

    node_selector = solver.node_selector
    components = _connected_components(A_sym)
    ordering = Vector{Int}(undef, n)
    k = 1

    for component in components
        submatrix = A_sym[component, component]
        component_ordering = _connected_cuthill_mckee_ordering(submatrix, node_selector)
        ordering[k:(k += length(component) - 1)] = component[component_ordering]
    end

    return ordering
end

"""
    _connected_cuthill_mckee_ordering(A, node_selector) -> Vector{Int}

TODO: Write here
"""
function _connected_cuthill_mckee_ordering(A::AbstractMatrix{Bool}, node_selector::Function)
    n = size(A, 1)
    ordering = Vector{Int}(undef, n)

    start = node_selector(A)
    degrees = vec(sum(A; dims=1))
    visited = Set(start)
    queue = Queue{Int}()
    enqueue!(queue, start)

    for i in 1:n
        parent = dequeue!(queue)
        ordering[i] = parent

        neighbors = findall(A[:, parent])
        unvisited = filter(!in(visited), neighbors)
        sort!(unvisited; by=node -> degrees[node])

        union!(visited, unvisited)
        foreach(neighbor -> enqueue!(queue, neighbor), unvisited)
    end

    return ordering
end

"""
    _connected_components(A) -> Vector{Vector{Int}}

TODO: Write here
"""
function _connected_components(A::AbstractMatrix{Bool})
    n = size(A, 1)
    visited = falses(n)
    queue = Queue{Int}()
    components = Vector{Int}[]

    for i in 1:n
        if !visited[i]
            visited[i] = true
            enqueue!(queue, i)
            component = Int[]

            while !isempty(queue)
                u = dequeue!(queue)
                push!(component, u)

                for v in 1:n
                    if A[u, v] && !visited[v]
                        visited[v] = true
                        enqueue!(queue, v)
                    end
                end
            end

            push!(components, component)
        end
    end

    return components
end
