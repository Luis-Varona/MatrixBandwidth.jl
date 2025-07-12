# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    DelCorsoManzini <: AbstractDecider <: AbstractAlgorithm

The *Del Corso–Manzini recognition algorithm* is a method for determining, given some fixed
``k ∈ ℕ``, whether a structurally symmetric matrix ``A`` has bandwidth at most ``k`` up to
symmetric permutation. The algorithm performs a depth-first search of all partial orderings
of the rows and columns of ``A``, adding indices one at a time. Partial orderings are pruned
not only by ensuring that adjacent pairs of currently placed indices are within ``k`` of
each other but also by tracking the latest positions at which the remaining indices can be
placed [DCM99](@cite).

As noted above, the Del Corso–Manzini algorithm requires structurally symmetric input (that
is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for ``1 ≤ i, j ≤ n``).

# Performance
Given an ``n×n`` input matrix ``A`` and threshold bandwidth ``k``, the Del Corso–Manzini
algorithm runs in ``O(nk ⋅ n!)`` time:
- We perform a depth-first search of ``O(n!)`` partial orderings.
- Checking plausibility of each partial ordering takes ``O(nk)`` time.
- Therefore, the overall time complexity is ``O(nk ⋅ n!)``.

Of course, this is but an upper bound on the time complexity of Del Corso–Manzini, achieved
only in the most pathological of cases. In practice, efficient pruning techniques and
compatibility checks result in approximately exponential growth in time complexity with
respect to ``n``.

Based on experimental results, the algorithm is feasible for ``n×n`` matrices up to
``n ≈ 100`` or so.

# Examples
[TODO: Write here]

# Notes
For readers of the original paper, what we call the Del Corso–Manzini recognition algorithm
here is essentially a wrapper around the underlying `AddNode` subroutine in what
[DCM99; p. 191](@cite) terms the "MB-ID algorithm" for bandwidth minimization (not mere
recognition). MB-ID (which we also implement in
[`MatrixBandwidth.Minimization.Exact.DelCorsoManzini`](@ref)) calls this recognition
procedure with incrementing values of ``k`` until a bandwidth-``k`` ordering is found, with
``k`` initialized to some lower bound on the minimum bandwidth of ``A`` up to symmetric
permutation.

[DCM99; p. 193](@cite) also describes an "MB-PS algorithm" for bandwidth minimization,
which we implement in [`MatrixBandwidth.Minimization.Exact.DelCorsoManziniWithPS`](@ref).
Similarly, the underlying recognition subroutine for MB-PS is implemented in
[`DelCorsoManziniWithPS`](@ref).
"""
struct DelCorsoManzini <: AbstractDecider end

Base.summary(::DelCorsoManzini) = "Del Corso–Manzini"

_requires_symmetry(::DelCorsoManzini) = true

"""
    DelCorsoManziniWithPS{D} <: AbstractDecider <: AbstractAlgorithm

The *Del Corso–Manzini recognition algorithm with perimeter search* is a method for
determining, given some fixed ``k ∈ ℕ``, whether a structurally symmetric matrix ``A`` has
bandwidth at most ``k`` up to symmetric permutation. The base Del Corso–Manzini algorithm
performs a depth-first search of all partial orderings of the rows and columns of ``A``,
adding indices one at a time. Partial orderings are pruned not only by ensuring that
adjacent pairs of currently placed indices are within ``k`` of each other but also by
tracking the latest positions at which the remaining indices can be placed [DCM99](@cite).

The incorporation of perimeter search to this approach entails precomputing a "perimeter" of
``d``-permutations of row indices of ``A``, where ``d`` is a positive integer passed as a
parameter to the solver. Each permutation represents a way to select the last ``d`` entries
of the ordering, and as the construction of the partial ordering progresses, potential
endings are pruned to exclude those incompatible with already placed indices. In addition to
pruning a potential ending if it contains indices already placed, compatibility is also
checked via precomputed time stamps indicating, for each potential ending, a loose lower
bound on the earliest position at which any given index can be placed should said ending be
selected.

As noted above, the Del Corso–Manzini algorithm with perimeter search requires structurally
symmetric input (that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is
nonzero for ``1 ≤ i, j ≤ n``).

# Fields
[TODO: Write here]

# Constructors
[TODO: Write here]

# Performance
Given an ``n×n`` input matrix ``A``, perimeter search depth ``d``, and threshold bandwidth
``k``, the Del Corso–Manzini algorithm with perimeter search runs in ``O(max(nᵈ, nk) ⋅ n!)``
time:
- We perform a depth-first search of ``O(n!)`` partial orderings.
- Checking plausibility of each partial ordering takes ``O(nk)`` time, and checking
    compatibility with all size-``d`` LPOs takes ``O(nᵈ)`` time. Thus, the overall time
    complexity for each value of ``k`` is ``O((nᵈ + nk) ⋅ n!)``.
- Therefore, the overall time complexity is ``O(max(nᵈ, nk) ⋅ n!)``.

Of course, this is but an upper bound on the time complexity of Del Corso–Manzini with
perimeter search, achieved only in the most pathological of cases. In practice, efficient
pruning techniques and compatibility checks result in approximately exponential growth in
time complexity with respect to ``n``.

Based on experimental results, the algorithm is feasible for ``n×n`` matrices up to
``n ≈ 100`` or so.

# Examples
[TODO: Write here]

# Notes
For readers of the original paper, what we call the Del Corso–Manzini recognition algorithm
with perimeter search here is essentially a wrapper around the underlying `AddNode1` and
`Prune` subroutines in what [DCM99; p. 193](@cite) terms the "MB-PS algorithm" for bandwidth
minimization (not mere recognition). MB-PS (which we also implement in
[`MatrixBandwidth.Minimization.Exact.DelCorsoManziniWithPS`](@ref)) calls this recognition
procedure with incrementing values of ``k`` until a bandwidth-``k`` ordering is found, with
``k`` initialized to some lower bound on the minimum bandwidth of ``A`` up to symmetric
permutation.

[DCM99; p. 191](@cite) also describes an "MB-ID algorithm" for bandwidth minimization, which
we implement in [`MatrixBandwidth.Minimization.Exact.DelCorsoManzini`](@ref). Similarly, the
underlying recognition subroutine for MB-ID is implemented in [`DelCorsoManzini`](@ref).
"""
struct DelCorsoManziniWithPS{D<:Union{Nothing,Int}} <: AbstractDecider
    depth::D

    DelCorsoManziniWithPS() = new{Nothing}(nothing)

    function DelCorsoManziniWithPS(depth::Int)
        if depth <= 0
            throw(ArgumentError("Perimeter search depth must be positive, got $depth"))
        end

        return new{Int}(depth)
    end
end

Base.summary(::DelCorsoManziniWithPS) = "Del Corso–Manzini with perimeter search"

_requires_symmetry(::DelCorsoManziniWithPS) = true

# TODO: Add inline comments to all the code in the rest of the file below

function _dcm_ps_optimal_depth(A::AbstractMatrix{Bool})
    n = size(A, 1)

    #= `A` is symmetric and diagonal entries are irrelevant (and indeed all false in this
    case), so we need only count entries above the main diagonal. =#
    num_nonzero = 2count(A[i, j] for i in 1:(n - 1) for j in (i + 1):n)

    #= Compute a lower bound on the bandwidth up to symmetric permutation and hence an upper
    bound on the band density given a minimal bandwidth ordering. =#
    k_lb = bandwidth_lower_bound(A)
    band_density_ub = num_nonzero / (2n * k_lb - k_lb * (k_lb + 1))

    # TODO: Explain
    if n < 100
        max_depth = 1.6e-8n^5 - 3.52e-6n^4 + 2.33e-4n^3 - 5.5e-3n^2 + 0.14n + 1.0
    else
        max_depth = 1
    end

    #= The optimal depth formula is nonincreasing with respect to band density and
    performance degradation from overestimation is far more serious than missed
    optimizations from underestimation, so we use our upper bound on the band density. =#
    density_contribution = -4(log(2.25band_density_ub))
    optimal_depth = clamp(max_depth + density_contribution, 1, max_depth)

    return round(Int, optimal_depth)
end

function _bool_bandwidth_k_ordering(A::AbstractMatrix{Bool}, k::Int, ::DelCorsoManzini)
    n = size(A, 1)

    ordering_buf = Vector{Int}(undef, n)
    adj_lists = map(node -> findall(A[:, node]), 1:n)

    unselected = Set(1:n)
    adj_list = Set{Int}()
    perimeter = [(Int[], Int[])]  # Vector of tuples (partial ordering, time stamps)
    num_placed = 0
    ps_depth = 0

    return _dcm_add_node!(
        ordering_buf, A, k, adj_lists, unselected, adj_list, perimeter, num_placed, ps_depth
    )
end

function _bool_bandwidth_k_ordering(
    A::AbstractMatrix{Bool}, k::Int, ::DelCorsoManziniWithPS{Nothing}
)
    decider_with_depth = DelCorsoManziniWithPS(_dcm_ps_optimal_depth(A))
    return _bool_bandwidth_k_ordering(A, k, decider_with_depth)
end

function _bool_bandwidth_k_ordering(
    A::AbstractMatrix{Bool}, k::Int, decider::DelCorsoManziniWithPS{Int}
)
    n = size(A, 1)

    ordering_buf = Vector{Int}(undef, n)
    adj_lists = map(node -> findall(A[:, node]), 1:n)

    unselected = Set(1:n)
    adj_list = Set{Int}()
    num_placed = 0
    ps_depth = decider.depth
    perimeter = [
        (lpo, _dcm_lpo_time_stamps(lpo, A, k)) for
        node_subset in combinations(1:n, ps_depth) for lpo in permutations(node_subset)
    ]

    return _dcm_add_node!(
        ordering_buf, A, k, adj_lists, unselected, adj_list, perimeter, num_placed, ps_depth
    )
end

function _dcm_add_node!(
    ordering_buf::Vector{Int},
    A::AbstractMatrix{Bool},
    k::Int,
    adj_lists::Vector{Vector{Int}},
    unselected::Set{Int},
    adj_list::Set{Int},
    perimeter::Vector{Tuple{Vector{Int},Vector{Int}}},
    num_placed::Int,
    ps_depth::Int,
)
    if isempty(perimeter)
        return nothing
    end

    n = size(A, 1)

    #= The original Del Corso and Manzini (1999) paper makes two bold claims:
    - Once `num_placed` reaches `n - ps_depth`, only one LPO will remain in `perimeter`.
    - This LPO will, by construction, be compatible with the current partial ordering under
        the bandwidth-`k` constraint.

    On both counts, the paper is wrong. The second claim assumes perfect pruning of the
    perimeter at each step, but the paper itself admits that this is infeasible, ultimately
    espousing the simplified version implemented here. The first claim, on the other hand,
    is never justified by the paper under *any* conditions—after all, many partial orderings
    of length `n - ps_depth` may allow multiple permutations of the remaining `ps_depth`
    indices without violating the bandwidth constraint.

    Therefore, we refrain from automatically filling up the remaining positions and
    returning early when `num_placed` reaches `n - ps_depth` (as the paper suggests),
    instead opting to continue the search until `num_placed` reaches `n`. =#
    if num_placed == n
        return ordering_buf
    end

    ordering = nothing
    res = iterate(unselected)

    while (!isnothing(res) && isnothing(ordering))
        (candidate, state) = res

        #= Ensure that no placed node at least `k` positions before `candidate` is adjacent
        to `candidate` (otherwise, the bandwidth would be at least `k + 1`). =#
        if all(
            !A[candidate, far_node] for far_node in view(ordering_buf, 1:(num_placed - k))
        )
            unselected_new = setdiff(unselected, [candidate])

            if _dcm_no_ps_is_reversed(
                ordering_buf, unselected_new, num_placed, ps_depth, candidate
            )
                continue
            end

            perimeter_new = _dcm_pruned_perimeter(
                perimeter, candidate, num_placed, ps_depth, n
            )

            if !isempty(perimeter_new)
                adj_list_new = intersect(
                    union(adj_list, adj_lists[candidate]), unselected_new
                )

                if _dcm_is_compatible(ordering_buf, A, adj_list_new, k, num_placed)
                    ordering_buf[num_placed + 1] = candidate
                    ordering = _dcm_add_node!(
                        ordering_buf,
                        A,
                        k,
                        adj_lists,
                        unselected_new,
                        adj_list_new,
                        perimeter_new,
                        num_placed + 1,
                        ps_depth,
                    )
                end
            end
        end

        res = iterate(unselected, state)
    end

    return ordering
end

#= The ordering `i₁, i₂, … iₙ` induces the same bandwidth as `iₙ, iₙ₋₁, … i₁`, so without
loss of generality, we restrict our search to partial orderings such that `i₁` is less than
the maximum possible `iₙ`. =#
function _dcm_no_ps_is_reversed(
    ordering_buf::Vector{Int},
    unselected::Set{Int},
    num_placed::Int,
    ps_depth::Int,
    candidate::Int,
)
    # If perimeter search is enabled, this check is already handled in the pruning logic
    if ps_depth > 0
        return false
    end

    #= If `num_placed` is zero and we check `ordering_buf[1] < max_last_label`,
    we risk erroneous behavior if `ordering_buf[1]` is initialized to some large
    value (e.g., with `UndefInitializer()`). =#
    if num_placed == 0
        first_label = candidate
    else
        first_label = ordering_buf[1]
    end

    #= If `unselected` is empty, then `candidate` is already the last label, not some label
    placed at an intermediate position that we can ignore. =#
    if isempty(unselected)
        max_last_label = candidate
    else
        max_last_label = maximum(unselected)
    end

    # There is no need to check for equality, since every label is unique
    return max_last_label < first_label
end

function _dcm_lpo_time_stamps(lpo::Vector{Int}, A::AbstractMatrix{Bool}, k::Int)
    n = size(A, 1)

    time_stamps = zeros(Int, n)
    d = length(lpo)
    foreach(((i, node),) -> time_stamps[node] = n - d + i, enumerate(lpo))

    queue = Queue{Int}()
    foreach(node -> enqueue!(queue, node), lpo)

    while !isempty(queue)
        node = dequeue!(queue)
        unvisited = filter!(neighbor -> time_stamps[neighbor] == 0, findall(A[:, node]))
        time_stamps[unvisited] .= time_stamps[node] - k
        foreach(neighbor -> enqueue!(queue, neighbor), unvisited)
    end

    return time_stamps
end

function _dcm_pruned_perimeter(
    perimeter::Vector{Tuple{Vector{Int},Vector{Int}}},
    candidate::Int,
    num_placed::Int,
    search_depth::Int,
    n::Int,
)
    if search_depth == 0
        perimeter_new = perimeter
    else
        perimeter_new = copy(perimeter)

        #= The ordering `i₁, i₂, … iₙ` induces the same bandwidth as `iₙ, iₙ₋₁, … i₁`, so
        without loss of generality, we restrict our search to partial orderings such that
        `i₁` is less than `iₙ`. =#
        if num_placed == 0
            filter!(((lpo, _),) -> candidate < lpo[end], perimeter_new)
        end

        if num_placed >= n - search_depth
            #= We are already filling up the last `search_depth` positions, so we restrict
            the search to LPOs in which `candidate` is placed at the next position. =#
            rel = ==
        else
            #= We are still filling up the first `n - search_depth` positions, so we allow
            all LPOs that are compatible with `candidate` being placed anywhere at or before
            the next position. =#
            rel = <=
        end

        filter!(
            ((_, time_stamp),) -> rel(time_stamp[candidate], num_placed + 1), perimeter_new
        )
    end

    return perimeter_new
end

function _dcm_is_compatible(
    ordering::Vector{Int},
    A::AbstractMatrix{Bool},
    adj_list::Set{Int},
    k::Int,
    num_placed::Int,
)
    if length(adj_list) > k
        return false
    end

    l = length(adj_list)
    latest_positions = Vector{Int}(undef, l)

    for (i, neighbor) in enumerate(adj_list)
        latest_position = typemax(Int)

        for (j, placed_node) in enumerate(view(ordering, 1:num_placed))
            if A[neighbor, placed_node]
                latest_position = min(latest_position, k + j)
            end
        end

        latest_positions[i] = latest_position
    end

    sort!(latest_positions)
    constraints = (1:l) .+ num_placed

    return all(latest_positions .>= constraints)
end
