# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CapraraSalazarGonzalez <: AbstractDecider <: AbstractAlgorithm

The *Caprara–Salazar-González recognition algorithm* is a method for determining, given some
fixed ``k ∈ ℕ``, whether a structurally symmetric matrix ``A`` has a bandwidth at most ``k``
up to symmetric permutation. The algorithm performs a depth-first search of all partial
orderings of the rows and columns of ``A``, adding indices one at a time. Partial orderings
are pruned not only by ensuring that adjacent pairs of currently placed indices are within
``k`` of each other but also by employing a branch-and-bound framework with lower bounds on
bandwidth compatibility computed via integer linear programming relaxations. This search is
repeated with incrementing values of ``k`` until a bandwidth-``k`` ordering is found [CS05],
with ``k`` initialized to some lower bound on the minimum bandwidth of ``A`` up to symmetric
permutation.

As noted above, the Caprara–Salazar-González algorithm requires structurally symmetric input
(that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for
``1 ≤ i, j ≤ n``).

# Supertype Hierarchy
`CapraraSalazarGonzalez` <: [`AbstractDecider`](@ref) <: [`AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix and threshold bandwidth ``k``, the Caprara–Salazar-González
algorithm runs in ``O(n! ⋅ mn)`` time in the worst case:

- We perform a depth-first search of ``O(n!)`` partial orderings.
- At each search node, we compute first and last feasible positions for all free nodes in
    ``O(mn)`` time using the formulas from Propositions 11 and 14 of [CS05].

Of course, this is all but an upper bound on the time complexity of
Caprara–Salazar-González, achieved only in the most pathological of cases. In practice,
efficient pruning techniques and compatibility checks result in approximately exponential
growth in time complexity with respect to ``n``.

# Examples
```@repl
using Random, SparseArrays
Random.seed!(17);
(n, p) = (10, 0.17);
A = sprand(n, n, p);
A = A + A' # Ensure structural symmetry;
(k_false, k_true) = (3, 6);
has_bandwidth_k_ordering(A, k_false, Recognition.CapraraSalazarGonzalez())
has_bandwidth_k_ordering(A, k_true, Recognition.CapraraSalazarGonzalez())
```

# Notes
For readers of the original paper, what we call the Caprara–Salazar-González algorithm here
is designated the `LAYOUT_LEFT_TO_RIGHT` algorithm in [CS05]. The paper also describes a
`LAYOUT_BOTH_WAYS` algorithm that performs a bidirectional search by adding indices to both
the left and right ends of the current partial ordering. However, this version is
considerably more complex to implement, and we ran into problems enforcing ILP constraints
on node pairs added to opposite ends of the ordering. In any case, computational results
demonstrate that neither `LAYOUT_LEFT_TO_RIGHT` nor `LAYOUT_BOTH_WAYS` is consistently
faster, and the paper states that there is no known heuristic for determining which version
will be more performant for a given input [CS05, pp. 368--69]. Therefore, we opt to
implement only `LAYOUT_LEFT_TO_RIGHT` as a matter of practicality, although future
developers may wish to extend the interface with `LAYOUT_BOTH_WAYS` as well.

Do also note that this algorithm is not the main `LAYOUT_LEFT_TO_RIGHT` procedure described
in the original paper, which actually never explicitly tackles matrix bandwidth recognition
[CS05]. However, the `LAYOUT_LEFT_TO_RIGHT` algorithm presented therein for bandwidth
*minimization* does repeatedly call a recognition subroutine—this is precisely the logic we
implement here. (We do, however, also implement said minimization algorithm in
[`MatrixBandwidth.Minimization.Exact.CapraraSalazarGonzalez`](@ref).)

# References
- [CS05](@cite): A. Caprara and J.-J. Salazar-González. *Laying Out Sparse Graphs with
    Provably Minimum Bandwidth*. INFORMS Journal on Computing **17**, 356–73 (2005).
    https://doi.org/10.1287/ijoc.1040.0083.
"""
struct CapraraSalazarGonzalez <: AbstractDecider end

push!(MatrixBandwidth.ALGORITHMS[:Recognition], CapraraSalazarGonzalez)

Base.summary(::CapraraSalazarGonzalez) = "Caprara–Salazar-González"

MatrixBandwidth._requires_structural_symmetry(::CapraraSalazarGonzalez) = true

function _has_bandwidth_k_ordering_impl(
    A::AbstractMatrix{Bool}, k::Integer, ::CapraraSalazarGonzalez
)
    n = size(A, 1)

    ordering_buf = Vector{Int}(undef, n)
    dist_matrix = floyd_warshall_shortest_paths(A)
    fixed = Int[]
    unselected = Set(1:n)

    return _csg_layout_left_to_right!(ordering_buf, A, k, dist_matrix, fixed, unselected)
end

function _csg_layout_left_to_right!(
    ordering_buf::Vector{Int},
    A::AbstractMatrix{Bool},
    k::Integer,
    dist_matrix::Matrix{Float64},
    fixed::Vector{Int},
    unselected::Set{Int},
)
    n = size(A, 1)

    if isempty(unselected)
        foreach(((i, node),) -> ordering_buf[i] = node, enumerate(fixed))
        return ordering_buf
    end

    earliest_positions, latest_positions = _csg_compute_positions(
        A, k, dist_matrix, fixed, unselected
    )

    if !_csg_feasible_positions(earliest_positions, latest_positions, fixed, unselected, n)
        return nothing
    end

    candidates = filter(v -> earliest_positions[v] == length(fixed) + 1, unselected)
    res = iterate(candidates)
    ordering = nothing

    while (!isnothing(res) && isnothing(ordering))
        candidate, state = res
        fixed_new = vcat(fixed, candidate)
        unselected_new = setdiff(unselected, candidate)

        ordering = _csg_layout_left_to_right!(
            ordering_buf, A, k, dist_matrix, fixed_new, unselected_new
        )

        res = iterate(candidates, state)
    end

    return ordering
end

function _csg_compute_positions(
    A::AbstractMatrix{Bool},
    k::Integer,
    dist_matrix::Matrix{Float64},
    fixed::Vector{Int},
    unselected::Set{Int},
)
    n = size(A, 1)
    fixed_set = Set(fixed) # For multiple `O(1)` membership checks below
    dists = Dict{Int,Float64}()

    for node in 1:n
        if node in fixed_set
            dists[node] = 0.0
        elseif isempty(fixed)
            dists[node] = Inf
        else
            dists[node] = minimum(dist_matrix[u, node] for u in fixed)
        end
    end

    earliest_positions = Dict{Int,Int}()
    latest_positions = Dict(u => i for (i, u) in enumerate(fixed))

    #= We process in order of increasing distance so closer neighbors (used to compute
    latest positions) have their latest positions computed first. =#
    for v in sort!(collect(unselected); by=v -> dists[v])
        dist_v = dists[v]
        lower_bound = length(fixed) + 1

        for (i, u) in enumerate(fixed)
            dist_uv = dist_matrix[u, v]

            # v must be within `k ⋅ dist(u, v)` positions of each fixed node `u`
            if isfinite(dist_uv)
                lower_bound = max(lower_bound, i - Int(dist_uv) * k)
            end
        end

        earliest_positions[v] = lower_bound
        neighbors = Int[]

        if isfinite(dist_v)
            for u in 1:n
                if v != u && A[v, u] && dists[u] == dist_v - 1
                    push!(neighbors, u)
                end
            end
        end

        if isinf(dist_v) || isempty(neighbors)
            latest_positions[v] = n
        else
            sorted_neighbors = sort(neighbors; by=u -> latest_positions[u], rev=true)
            #= `running_pos` tracks the last slot where all neighbors can occupy distinct
            positions. =#
            running_pos = latest_positions[sorted_neighbors[1]]

            for u in Iterators.drop(sorted_neighbors, 1)
                running_pos = min(latest_positions[u], running_pos - 1)
            end

            # `v` must be within `k` positions of its closest neighbors
            latest_positions[v] = k + running_pos
        end
    end

    return (earliest_positions, latest_positions)
end

function _csg_feasible_positions(
    earliest_positions::Dict{Int,Int},
    latest_positions::Dict{Int,Int},
    fixed::Vector{Int},
    unselected::Set{Int},
    n::Int,
)
    remaining = copy(unselected)
    i = length(fixed) + 1
    feasible = true

    while (i <= n && feasible)
        eligible = Iterators.filter(u -> earliest_positions[u] <= i, remaining)

        if isempty(eligible)
            feasible = false
        else
            # Greedily assign the most constrained eligible node to this position
            v = argmin(u -> latest_positions[u], eligible)

            if latest_positions[v] < i
                feasible = false
            else
                delete!(remaining, v)
            end
        end

        i += 1
    end

    return feasible
end
