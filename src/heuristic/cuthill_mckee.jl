# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CuthillMcKee <: HeuristicSolver <: AbstractSolver

The *Cuthill–McKee algorithm* is a heuristic method for minimizing the bandwidth of a
symmetric matrix ``A``. It considers the graph ``G(A)`` whose adjacency matrix is ``A``
(ignoring self-loops) and performs a breadth-first search of each connected component of
``G(A)``, starting from a low-degree node then visiting its neighbors in order of increasing
degree. Particularly effective when ``A`` is sparse, this heuristic typically produces an
ordering which induces a matrix bandwidth either equal to or very close to the true minimum
[CM69; pp. 157--58](@cite).

We also extend the algorithm to work more generally when ``A`` is not symmetric by applying
it to ``A + Aᵀ`` instead, as suggested in [RS06; p. 808](@cite). This approach still tends
to produce a fairly good ordering, but it is not guaranteed to be as optimal as directly
applying Cuthill–McKee to a symmetric input.

# Fields
- `node_selector::Function`: a function that selects a node from some connected component of
    the input matrix from which to start the breadth-first search. If no custom heuristic is
    specified, this field defaults to [`pseudo_peripheral_node`](@ref), which picks a node
    "farthest" from the others in the component (not necessarily the lowest-degree node).

# Examples
Cuthill–McKee finds an optimal ordering for an asymmetric ``35×35`` matrix with bandwidth
``3`` whose rows and columns have been shuffled:
```jldoctest
julia> using Random

julia> Random.seed!(13);

julia> (n, k) = (35, 3);

julia> perm = randperm(n);

julia> A = MatrixBandwidth.random_sparse_banded_matrix(n, k);

julia> A_shuffled = A[perm, perm];

julia> res = minimize_bandwidth(A_shuffled, CuthillMcKee());

julia> iszero.(A_shuffled) == iszero.(A_shuffled') # Works even for asymmetric matrices
false

julia> bandwidth_unpermuted(A)
3

julia> bandwidth_unpermuted(A_shuffled)
33

julia> res.bandwidth # The true minimum bandwidth
3
```

Cuthill–McKee finds a near-optimal ordering for an asymmetric ``200×200`` matrix with
bandwidth ``10`` whose rows and columns have been shuffled:
```jldoctest
julia> using Random

julia> Random.seed!(37452);

julia> (n, k) = (200, 10);

julia> perm = randperm(n);

julia> A = MatrixBandwidth.random_sparse_banded_matrix(n, k);

julia> A_shuffled = A[perm, perm];

julia> res = minimize_bandwidth(A_shuffled, CuthillMcKee());

julia> iszero.(A_shuffled) == iszero.(A_shuffled') # Works even for asymmetric matrices
false

julia> bandwidth_unpermuted(A)
10

julia> bandwidth_unpermuted(A_shuffled)
194

julia> res.bandwidth # Close to the true minimum
14
```

# Notes
It was found in [Geo71; pp. 114--15](@cite) that reversing the ordering produced by
Cuthill–McKee tends to produce a better ordering; this variant is known as the *reverse
Cuthill–McKee algorithm*. See [`ReverseCuthillMcKee`](@ref) and the associated method of
`_bool_minimal_band_ordering` for our implementation of this algorithm.

Note also that the `node_selector` field must be of the form
`(A::AbstractMatrix{Bool}) -> Integer` (i.e., it must take in an boolean matrix and return
an integer). If this is not the case, an `ArgumentError` is thrown upon construction.

See also the documentation for supertypes [`HeuristicSolver`](@ref) and
[`AbstractSolver`](@ref).
"""
struct CuthillMcKee <: HeuristicSolver
    node_selector::Function

    function CuthillMcKee(node_selector::Function=DEFAULT_SELECTOR)
        _assert_valid_node_selector(node_selector)
        return new(node_selector)
    end
end

Base.summary(::CuthillMcKee) = "Cuthill–McKee algorithm"

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::CuthillMcKee)
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

# Cuthill–McKee performs a breadth-first search on each connected component independently
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

# Find the indices of all connected components in an adjacency matrix
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
