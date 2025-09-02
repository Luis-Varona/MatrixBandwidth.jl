# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SaxeGurariSudborough <: AbstractDecider <: AbstractAlgorithm

[TODO: Write here]

# Supertype Hierarchy
`SaxeGurariSudborough` <: [`AbstractDecider`](@ref) <: [`AbstractAlgorithm`](@ref)

[TODO: Write here]

# Performance
Given an ``n×n`` input matrix ``A`` and threshold bandwidth ``k``, the
Saxe–Gurari–Sudborough algorithm runs in ``O(nᵏ)`` time [GS84, p. 531]. This is an
improvement upon the original ``O(nᵏ⁺¹)`` Saxe algorithm [Sax80, p. 363].

Of course, when ``k < 3``, then the initial ``O(n³)`` bandwidth lower bound computation
(specifically, [`bandwidth_lower_bound`](@ref)) performed in all
[`has_bandwidth_k_ordering`](@ref) calls dominates the overall complexity (although the
constant scaling factor of that subroutine is generally much smaller than that of the
algorithm proper).

# Examples
[TODO: Write here]

# Notes
Whereas most bandwidth recognition algorithms are technically factorial-time (with respect
to ``n``) in the worst case but practically always approximate exponential time complexity
in real life (see: [`DelCorsoManzini`](@ref)), the ``O(nᵏ)`` upper bound on
Saxe–Gurari–Sudborough is typically a good representation of actual performance in most
cases. Indeed, these other types of algorithms tend to outperform Saxe–Gurari–Sudborough
for larger ``k``, given that their aggressive pruning strategies keep their effective search
space very small in practice.

# References
- [GS84](@cite): E. M. Gurari and I. H. Sudborough. *Improved dynamic programming algorithms
    for bandwidth minimization and the MinCut Linear Arrangement problem*. Journal of
    Algorithms **5**, 531–46 (1984). https://doi.org/10.1016/0196-6774(84)90006-3.
- [Sax80](@cite): J. B. Saxe. *Dynamic-Programming Algorithms for Recognizing
    Small-Bandwidth Graphs in Polynomial Time*. SIAM Journal on Algebraic and Discrete
    Methods **1**, 363–69 (1980). https://doi.org/10.1137/0601042.
"""
struct SaxeGurariSudborough <: AbstractDecider end

push!(ALGORITHMS[:Recognition], SaxeGurariSudborough)

Base.summary(::SaxeGurariSudborough) = "Saxe–Gurari–Sudborough"

_requires_structural_symmetry(::SaxeGurariSudborough) = true

function _bool_bandwidth_k_ordering(
    A::AbstractMatrix{Bool}, k::Integer, ::SaxeGurariSudborough
)
    components = _connected_components(A)
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

# TODO: After this point, more comprehensive inline comments are required.

function _sgs_connected_ordering(A::AbstractMatrix{Bool}, k::Integer)
    n = size(A, 1)
    adj_lists = map(node -> findall(view(A, :, node)), 1:n)

    KeyType = Tuple{Tuple{Vararg{Int}},Tuple{Vararg{Tuple{Int,Int}}}}
    parent = Dict{KeyType,Tuple{Union{Nothing,KeyType},Union{Nothing,Int}}}()
    nums_placed = Dict{KeyType,Int}()
    queue = Queue{KeyType}()
    visited = Set{KeyType}()

    empty_key = ((), ())
    parent[empty_key] = (nothing, nothing)
    nums_placed[empty_key] = 0
    enqueue!(queue, empty_key)
    push!(visited, empty_key)

    while !isempty(queue)
        key = dequeue!(queue)
        num_placed = nums_placed[key]
        region = collect(Int, key[1])
        dangling = Set{Tuple{Int,Int}}(key[2])

        if length(region) < k
            candidates = _sgs_unassigned(region, dangling, adj_lists)
        else
            u = region[1]
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
                key_new = _sgs_key(region_new, dangling_new)

                if _sgs_layout_is_plausible(region_new, dangling_new, k)
                    num_placed_new = num_placed + 1

                    # if !(key_new in visited)
                    #     parent[key_new] = (key, v)
                    #     nums_placed[key_new] = num_placed_new
                    #     enqueue!(queue, key_new)
                    #     push!(visited, key_new)
                    # end

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

                    # if num_placed_new == n
                    #     ordering = Vector{Int}(undef, n)
                    #     (key_new, v) = parent[key_new]

                    #     while !isnothing(key_new)
                    #         ordering[(num_placed_new -= 1) + 1] = v
                    #         (key_new, v) = parent[key_new]
                    #     end

                    #     return ordering
                    # end

                    if !(key_new in visited)
                        parent[key_new] = (key, v)
                        nums_placed[key_new] = num_placed_new
                        enqueue!(queue, key_new)
                        push!(visited, key_new)
                    end
                end
            end
        end
    end

    return nothing
end

function _sgs_unassigned(
    region::Vector{Int}, dangling::Set{Tuple{Int,Int}}, adj_lists::Vector{Vector{Int}}
)
    if isempty(dangling)
        return Set(eachindex(adj_lists))
    end

    unassigned = Set{Int}()
    processed = Set{Tuple{Int,Int}}(dangling)
    queue = Queue{Tuple{Int,Int}}()
    foreach(edge -> enqueue!(queue, edge), dangling)
    visited = Set(region)

    @inline function _update_state!(node::Int)
        push!(unassigned, node)
        push!(visited, node)

        for neighbor in adj_lists[node]
            edge = _pot_edge(node, neighbor)

            if !(edge in processed)
                push!(processed, edge)
                enqueue!(queue, edge)
            end
        end
    end

    while !isempty(queue)
        (u, v) = dequeue!(queue)

        if !(u in visited)
            _update_state!(u)
        elseif !(v in visited)
            _update_state!(v)
        end
    end

    return unassigned
end

function _sgs_dangling_new(
    region::Vector{Int}, dangling::Set{Tuple{Int,Int}}, v::Int, adj_list::Vector{Int}
)
    neighbors_in_region = Iterators.filter(u -> _pot_edge(u, v) in dangling, region)
    dangling_new = setdiff(
        dangling, Iterators.map(u -> _pot_edge(u, v), neighbors_in_region)
    )
    additional_edges = Iterators.map(
        u -> _pot_edge(u, v), Iterators.filter(u -> u != v && !(u in region), adj_list)
    )
    return union!(dangling_new, additional_edges)
end

function _sgs_region_new(
    region_extended::Vector{Int}, dangling_new::Set{Tuple{Int,Int}}, k::Integer
)
    first_idx_in_region = Dict{Int,Int}()
    foreach(
        ((i, v),) -> first_idx_in_region[v] = get(first_idx_in_region, v, i),
        enumerate(region_extended),
    )
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

@inline function _sgs_key(region::Vector{Int}, dangling::Set{Tuple{Int,Int}})
    return (Tuple(region), Tuple(sort!(collect(dangling))))
end

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

@inline function _pot_edge(u::Int, v::Int)
    return (min(u, v), max(u, v))
end
