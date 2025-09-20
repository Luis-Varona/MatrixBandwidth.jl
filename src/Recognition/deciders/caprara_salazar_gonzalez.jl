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
algorithm runs in ``O(n! ⋅ Tᵢₗₚ(n, n²))`` time, where ``Tᵢₗₚ(n, m)`` is the time taken to
solve an integer linear programming (ILP) problem with ``O(n)`` variables and ``O(m)``
constraints:

- We perform a depth-first search of ``O(n!)`` partial orderings.
- At each search node, we solve ILP relaxations with ``n`` variables and ``O(n²)``
    constraints (given by the number of nonzero entries in the computed distance matrix),
    taking ``Tᵢₗₚ(n, n²)`` time. (This dominates the ``O(n²)`` auxiliary computations needed
    to set up the ILP.)
- Therefore, the overall time complexity is ``O(n! ⋅ Tᵢₗₚ(n, n²))``.

Note that ``Tᵢₗₚ(n, n²)`` has worst-case complexity ``O(2ⁿ)``, although this ultimately
depends on the ILP solver used. (Here, we use the HiGHS solver from the `HiGHS.jl` package.)

Of course, this is all but an upper bound on the time complexity of
Caprara–Salazar-González, achieved only in the most pathological of cases. In practice,
efficient pruning techniques and compatibility checks result in approximately exponential
growth in time complexity with respect to ``n``.

# Examples
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(17);

julia> (n, p) = (10, 0.17);

julia> A = sprand(n, n, p);

julia> A = A + A' # Ensure structural symmetry;

julia> (k_false, k_true) = (3, 6);

julia> has_bandwidth_k_ordering(A, k_false, Recognition.CapraraSalazarGonzalez())
Results of Bandwidth Recognition Algorithm
 * Algorithm: Caprara–Salazar-González
 * Bandwidth Threshold k: 3
 * Has Bandwidth ≤ k Ordering: false
 * Original Bandwidth: 9
 * Matrix Size: 10×10

julia> has_bandwidth_k_ordering(A, k_true, Recognition.CapraraSalazarGonzalez())
Results of Bandwidth Recognition Algorithm
 * Algorithm: Caprara–Salazar-González
 * Bandwidth Threshold k: 6
 * Has Bandwidth ≤ k Ordering: true
 * Original Bandwidth: 9
 * Matrix Size: 10×10
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

A final implementation detail worth noting is that we use HiGHS as our solver; it is one of
the fastest open-source solvers available for mixed-integer linear programming.

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

    earliest_positions = _csg_solve_relaxation(A, k, dist_matrix, fixed, unselected)

    if isnothing(earliest_positions) ||
        !_csg_feasible_positions(earliest_positions, fixed, unselected, n)
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

#= Use integer-linear programming to determine the earliest feasible position of each node
not yet placed in the ordering. =#
function _csg_solve_relaxation(
    A::AbstractMatrix{Bool},
    k::Integer,
    dist_matrix::Matrix{Float64},
    fixed::Vector{Int},
    unselected::Set{Int},
)
    earliest_positions = Dict{Int,Int}()
    dists = Dict{Int,Float64}()

    for node in unselected
        if isempty(fixed)
            dists[node] = Inf
        else
            dists[node] = minimum(dist_matrix[fixed, node])
        end
    end

    computed_earliest = Dict{Int,Int}()

    for node in sort!(collect(unselected); by=v -> dists[v])
        earliest_val = _csg_solve_inner_ilp(
            A, k, dist_matrix, fixed, computed_earliest, node
        )

        if isnothing(earliest_val)
            return nothing
        end

        earliest_positions[node] = earliest_val
        computed_earliest[node] = earliest_val
    end

    return earliest_positions
end

function _csg_solve_inner_ilp(
    A::AbstractMatrix{Bool},
    k::Integer,
    dist_matrix::Matrix{Float64},
    fixed::Vector{Int},
    computed_earliest::Dict{Int,Int},
    v::Int,
)
    n = size(A, 1)

    #= `HiGHS` is one of the fastest open-source solvers available for mixed-integer linear
    programming. =#
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, pos_v, Int)
    @objective(model, Min, pos_v)

    #= The next placed node must be placed somewhere after all currently fixed nodes, but
    (naturally) before the end of the ordering. =#
    @constraint(model, pos_v >= length(fixed) + 1)
    @constraint(model, pos_v <= n)

    for (i, u) in enumerate(fixed)
        dist_temp = dist_matrix[u, v]

        if isfinite(dist_temp)
            dist = Int(dist_temp)
            #= The difference in indices of `u` and `v` must be at most `k` times the
            shortest-path distance between them. If `u` and `v` have distance 1, this
            simplifies to the standard bandwidth constraint of adjacent nodes being at most
            `k` indices apart. =#
            @constraint(model, pos_v <= k * dist + i)
            @constraint(model, pos_v >= i - k * dist)
        end
    end

    for (u, pos_u) in computed_earliest
        dist_temp = dist_matrix[u, v]

        if isfinite(dist_temp)
            dist = Int(dist_temp)
            #= The difference in indices of `u` and `v` must be at most `k` times the
            shortest-path distance between them. If `u` and `v` have distance 1, this
            simplifies to the standard bandwidth constraint of adjacent nodes being at most
            `k` indices apart. =#
            @constraint(model, pos_v >= pos_u - k * dist)
            @constraint(model, pos_v <= pos_u + k * dist)
        end
    end

    optimize!(model)

    if termination_status(model) == OPTIMAL
        solution = Int(value(pos_v))
    else
        solution = nothing
    end

    return solution
end

function _csg_feasible_positions(
    earliest_positions::Dict{Int,Int}, fixed::Vector{Int}, unselected::Set{Int}, n::Int
)
    remaining = copy(unselected)

    for i in (length(fixed) + 1):n
        eligible_nodes = Iterators.filter(v -> earliest_positions[v] <= i, remaining)

        if isempty(eligible_nodes)
            return false
        end

        delete!(remaining, first(eligible_nodes))
    end

    return true
end
