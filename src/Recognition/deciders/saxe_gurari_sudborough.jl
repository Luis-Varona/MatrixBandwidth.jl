# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SaxeGurariSudborough <: AbstractDecider <: AbstractAlgorithm

The *Saxe–Gurari–Sudborough recognition algorithm* is a method for determining, given some
fixed ``k ∈ ℕ``, whether a structurally symmetric matrix ``A`` has bandwidth at most ``k``
up to symmetric permutation. The algorithm employs dynamic programming to search over
equivalence classes of partial orderings, where two partial orderings of length ``l`` are
equivalent if they share the same *active region*. (The active region of a partial ordering
is defined as the sequence of the last ``min(k, l)`` vertices in the ordering taken together
with all *dangling edges*—edges with one endpoint in the ordering and the other endpoint not
yet in the ordering.) It extends these partial layouts one vertex at a time in a
breadth-first manner, pruning implausible classes that violate bandwidth-``k`` constraints
such as degree bounds on active vertices and excessive numbers of dangling edges [GS84].

As noted above, the Saxe–Gurari–Sudborough algorithm requires structurally symmetric input
(that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for
``1 ≤ i, j ≤ n``).

# Supertype Hierarchy
`SaxeGurariSudborough` <: [`AbstractDecider`](@ref) <: [`AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix and threshold bandwidth ``k``, the Saxe–Gurari–Sudborough
algorithm runs in ``O(nᵏ)`` time [GS84, p. 531]. This is an improvement upon the original
``O(nᵏ⁺¹)`` Saxe algorithm [Sax80, p. 363]. (Of course, when ``k < 3``, then the initial
``O(n³)`` bandwidth lower bound computation performed in all
[`has_bandwidth_k_ordering`](@ref) calls dominates the overall complexity, although the
constant scaling factor of that subroutine is generally much smaller than that of the
algorithm proper).

Whereas most bandwidth recognition algorithms are technically factorial-time (with respect
to ``n``) in the worst case but practically always approximate exponential time complexity
in real life (see: [`DelCorsoManzini`](@ref)), the ``O(nᵏ)`` upper bound on
Saxe–Gurari–Sudborough is typically a good representation of actual performance in most
cases. Indeed, these other types of algorithms tend to outperform Saxe–Gurari–Sudborough
for larger ``k``, given that their aggressive pruning strategies keep their effective search
space very small in practice.

# Examples
```@repl
using Random, SparseArrays
Random.seed!(274);
(n, p) = (20, 0.08);
A = sprand(n, n, p);
A = A + A' # Ensure structural symmetry;
(k_false, k_true) = (3, 5);
has_bandwidth_k_ordering(A, k_false, Recognition.SaxeGurariSudborough())
has_bandwidth_k_ordering(A, k_true, Recognition.SaxeGurariSudborough())
```

# Notes
This general family of bandwidth recognition algorithms was conceived as a response to a
question posed by [GGJK78, p. 494]: is the "bandwidth ≤ k?" problem NP-complete for
arbitrary ``k``? [Sax80] answered this question in the negative by providing a ``O(nᵏ⁺¹)``
algorithm, constructively proving that the problem is class P. Later, [GS84] improved upon
this algorithm by reducing time complexity to ``O(nᵏ)``. Whereas the original Saxe algorithm
considers extensions of partial orderings with any remaining unplaced vertex (of which there
are ``O(n)`` at any point in the breadth-first search), the Gurari–Sudborough refinement
only considers extensions with vertices reachable by paths beginning with a dangling edge
that never again traverse a dangling edge [GS84, pp. 535–36].

# References
- [GGJK78](@cite): M. R. Garey, R. L. Graham, D. S. Johnson and D. E. Knuth. *Complexity
    Results for Bandwidth Minimization*. SIAM Journal on Applied Mathematics **34**, 477–95
    (1978). https://doi.org/10.1137/0134037.
- [GS84](@cite): E. M. Gurari and I. H. Sudborough. *Improved dynamic programming algorithms
    for bandwidth minimization and the MinCut Linear Arrangement problem*. Journal of
    Algorithms **5**, 531–46 (1984). https://doi.org/10.1016/0196-6774(84)90006-3.
- [Sax80](@cite): J. B. Saxe. *Dynamic-Programming Algorithms for Recognizing
    Small-Bandwidth Graphs in Polynomial Time*. SIAM Journal on Algebraic and Discrete
    Methods **1**, 363–69 (1980). https://doi.org/10.1137/0601042.
"""
struct SaxeGurariSudborough <: AbstractDecider end

push!(MatrixBandwidth.ALGORITHMS[:Recognition], SaxeGurariSudborough)

Base.summary(::SaxeGurariSudborough) = "Saxe–Gurari–Sudborough"

MatrixBandwidth._requires_structural_symmetry(::SaxeGurariSudborough) = true

function _has_bandwidth_k_ordering_impl(
    A::AbstractMatrix{Bool}, k::Integer, ::SaxeGurariSudborough
)
    components = connected_components(A)
    #= Smaller components take less time to process, so we do them first to heuristically
    break earlier in some situtations and avoid wasting time on larger components. =#
    sort!(components; by=length)

    n = size(A, 1)
    ordering = Vector{Int}(undef, n)
    res = iterate(components)
    num_placed = 0

    while !isnothing(res)
        (component, state) = res
        submatrix = view(A, component, component)
        component_ordering = _sgs_connected_ordering(submatrix, k)

        if isnothing(component_ordering)
            res = nothing
        else
            ordering[(num_placed + 1):(num_placed += length(component))] .= component[component_ordering]
            res = iterate(components, state)
        end
    end

    if num_placed < n
        ordering = nothing
    end

    return ordering
end

# Saxe–Gurari–Sudborough searches each connected component independently
function _sgs_connected_ordering(A::AbstractMatrix{Bool}, k::Integer)
    n = size(A, 1)
    adj_lists = map(node -> findall(view(A, :, node)), 1:n)

    edge_idxs = Dict{Tuple{Int,Int},Int}()
    num_edges = 0

    for (u, neighbors) in enumerate(adj_lists), v in Iterators.filter(w -> w > u, neighbors)
        edge = (u, v)
        edge_idxs[edge] = (num_edges += 1)
    end

    KeyType = Tuple{Tuple{Vararg{Int}},Tuple{Vararg{Tuple{Int,Int}}}}
    parent = Dict{KeyType,Tuple{Union{Nothing,KeyType},Union{Nothing,Int}}}()
    nums_placed = Dict{KeyType,Int}()
    queue = Queue{KeyType}()
    visited = Set{KeyType}()

    empty_key = ((), ())
    parent[empty_key] = (nothing, nothing)
    nums_placed[empty_key] = 0
    push!(queue, empty_key)
    push!(visited, empty_key)

    nodes_visited = falses(n)
    edges_visited = falses(num_edges)

    while !isempty(queue)
        key = popfirst!(queue)
        num_placed = nums_placed[key]
        region = collect(Int, key[1])
        dangling = Set{Tuple{Int,Int}}(key[2])

        if length(region) < k
            #= The Gurari–Sudborough refinement of the original Saxe algorithm cuts time
            complexity down from `O(nᵏ⁺¹)` to `O(nᵏ)` by only considering extensions with
            vertices reachable by a path beginning with a dangling edge that never again
            traverses a dangling edge. =#
            candidates = _sgs_unassigned(
                region, dangling, adj_lists, nodes_visited, edges_visited, edge_idxs
            )
        else
            u = region[1]
            #= Indeed, there should be exactly one such edge; the lazy evaluation of
            `Iterators.filter` allows us to break once said edge is found. =#
            edge = first(Iterators.filter(edge -> u in edge, dangling))

            if u == edge[1]
                v = edge[2]
            else
                v = edge[1]
            end

            candidates = Set(v)
        end

        for v in candidates
            dangling_new = _sgs_dangling_new(region, dangling, v, adj_lists[v])
            region_extended = vcat(region, v)
            region_new = _sgs_region_new(region_extended, dangling_new, k)

            if !isnothing(region_new)
                if _sgs_layout_is_plausible(region_new, dangling_new, k)
                    key_new = _sgs_key(region_new, dangling_new)
                    num_placed_new = num_placed + 1

                    if num_placed_new == n
                        ordering = Vector{Int}(undef, n)
                        ordering[num_placed_new] = v
                        (key_new, v) = parent[key]

                        while !isnothing(key_new)
                            ordering[num_placed_new -= 1] = v
                            (key_new, v) = parent[key_new]
                        end

                        return ordering
                    end

                    if !(key_new in visited)
                        parent[key_new] = (key, v)
                        nums_placed[key_new] = num_placed_new
                        push!(queue, key_new)
                        push!(visited, key_new)
                    end
                end
            end
        end
    end

    return nothing
end

#= Return the set of vertices not yet in `region` reachable from `region` via paths
beginning with a dangling edge that never again traverse a dangling edge. (If `dangling` is
empty, then of course, all unassigned vertices are reachable.) =#
function _sgs_unassigned(
    region::Vector{Int},
    dangling::Set{Tuple{Int,Int}},
    adj_lists::Vector{Vector{Int}},
    nodes_visited::BitVector,
    edges_visited::BitVector,
    edge_idxs::Dict{Tuple{Int,Int},Int},
)
    if isempty(dangling)
        return Set(eachindex(adj_lists))
    end

    unassigned = Set{Int}()
    fill!(nodes_visited, false)
    fill!(edges_visited, false)
    foreach(edge -> edges_visited[edge_idxs[edge]] = true, dangling)

    queue = Queue{Tuple{Int,Int}}()
    foreach(edge -> push!(queue, edge), dangling)
    nodes_visited[region] .= true

    @inline function _update_state!(node::Int)
        push!(unassigned, node)
        nodes_visited[node] = true

        for neighbor in adj_lists[node]
            edge = _pot_edge(node, neighbor)
            edge_idx = edge_idxs[edge]

            if !edges_visited[edge_idx]
                edges_visited[edge_idx] = true
                push!(queue, edge)
            end
        end
    end

    while !isempty(queue)
        (u, v) = popfirst!(queue)

        if !nodes_visited[u]
            _update_state!(u)
        elseif !nodes_visited[v]
            _update_state!(v)
        end
    end

    return unassigned
end

# Update the set of dangling edges after adding `v` to `region`
function _sgs_dangling_new(
    region::Vector{Int}, dangling::Set{Tuple{Int,Int}}, v::Int, adj_list::Vector{Int}
)
    neighbors_in_region = Iterators.filter(u -> _pot_edge(u, v) in dangling, region)
    dangling_new = setdiff(
        dangling, Iterators.map(u -> _pot_edge(u, v), neighbors_in_region)
    )
    additional_edges = Iterators.map(
        u -> _pot_edge(u, v), Iterators.filter(!in(region), adj_list)
    )
    return union!(dangling_new, additional_edges)
end

# Update the active region after adding `v` to `region`
function _sgs_region_new(
    region_extended::Vector{Int}, dangling_new::Set{Tuple{Int,Int}}, k::Integer
)
    first_idx_in_region = Dict{Int,Int}()
    foreach(((i, v),) -> first_idx_in_region[v] = i, enumerate(region_extended))
    region_extended_set = Set(region_extended)
    active_new = Set(
        Iterators.filter(in(region_extended_set), Iterators.flatten(dangling_new))
    )

    if isempty(active_new)
        region_new = Int[]
    else
        extension_size = length(region_extended)
        pos_min = minimum(Iterators.map(v -> first_idx_in_region[v], active_new))

        if extension_size - pos_min + 1 > k
            region_new = nothing
        else
            region_new = region_extended[pos_min:extension_size]
        end
    end

    return region_new
end

#= Construct a hashable key for the current active region (taken together with the set of
dangling edges). =#
@inline function _sgs_key(region::Vector{Int}, dangling::Set{Tuple{Int,Int}})
    return (Tuple(region), Tuple(sort!(collect(dangling))))
end

#= Check whether the current layout is plausible given bandwidth-`k` constraints such as
degree bounds on active vertices and excessive numbers of dangling edges. =#
function _sgs_layout_is_plausible(
    region::Vector{Int}, dangling::Set{Tuple{Int,Int}}, k::Integer
)
    region_size = length(region)
    nums_dangling = Iterators.map(v -> count(edge -> v in edge, dangling), region)
    return all(
        ((num_dangling, idx),) -> num_dangling <= k - region_size + idx,
        zip(nums_dangling, 1:region_size),
    )
end

#= Represent the (potential) edge between `u` and `v` as a sorted tuple for hashing
purposes (as opposed to a Set). =#
@inline function _pot_edge(u::Int, v::Int)
    return (min(u, v), max(u, v))
end
