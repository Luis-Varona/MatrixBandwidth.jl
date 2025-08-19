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
placed [DM99].

As noted above, the Del Corso–Manzini algorithm requires structurally symmetric input (that
is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for ``1 ≤ i, j ≤ n``).

# Supertype Hierarchy
`DelCorsoManzini` <: [`AbstractDecider`](@ref) <: [`AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix ``A`` and threshold bandwidth ``k``, the Del Corso–Manzini
algorithm runs in ``O(n! ⋅ nk)`` time:
- We perform a depth-first search of ``O(n!)`` partial orderings.
- Checking plausibility of each partial ordering takes ``O(nk)`` time.
- Therefore, the overall time complexity is ``O(n! ⋅ nk)``.

Of course, this is but an upper bound on the time complexity of Del Corso–Manzini, achieved
only in the most pathological of cases. In practice, efficient pruning techniques and
compatibility checks result in approximately exponential growth in time complexity with
respect to ``n``.

Based on experimental results, the algorithm is feasible for ``n×n`` matrices up to
``n ≈ 100`` or so.

# Examples
We demonstrate both an affirmative and a negative result for the Del Corso–Manzini
recognition algorithm on a random ``40×40`` matrix:
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(7878);

julia> (n, p) = (40, 0.1);

julia> A = sprand(n, n, p);

julia> A = A + A' # Ensure structural symmetry;

julia> (k_false, k_true) = (13, 26);

julia> has_bandwidth_k_ordering(A, k_false, Recognition.DelCorsoManzini())
Results of Bandwidth Recognition Algorithm
 * Algorithm: Del Corso–Manzini
 * Bandwidth Threshold k: 13
 * Has Bandwidth ≤ k Ordering: false
 * Original Bandwidth: 34
 * Matrix Size: 40×40

julia> has_bandwidth_k_ordering(A, k_true, Recognition.DelCorsoManzini())
Results of Bandwidth Recognition Algorithm
 * Algorithm: Del Corso–Manzini
 * Bandwidth Threshold k: 26
 * Has Bandwidth ≤ k Ordering: true
 * Original Bandwidth: 34
 * Matrix Size: 40×40
```

# Notes
For readers of the original paper, what we call the Del Corso–Manzini recognition algorithm
here is essentially a wrapper around the underlying `AddNode` subroutine in what
[DM99, p. 191] term the "MB-ID algorithm" for bandwidth minimization (not mere recognition).
MB-ID (which we also implement in
[`MatrixBandwidth.Minimization.Exact.DelCorsoManzini`](@ref)) calls this recognition
procedure with incrementing values of ``k`` until a bandwidth-``k`` ordering is found, with
``k`` initialized to some lower bound on the minimum bandwidth of ``A`` up to symmetric
permutation.

[DM99, p. 193] also describe an "MB-PS algorithm" for bandwidth minimization, which we
implement in [`MatrixBandwidth.Minimization.Exact.DelCorsoManziniWithPS`](@ref). Similarly,
the underlying recognition subroutine for MB-PS is implemented in
[`DelCorsoManziniWithPS`](@ref).

# References
- [DM99](@cite): G. M. Del Corso and G. Manzini. *Finding Exact Solutions to the Bandwidth
    Minimization Problem*. Computing **62**, 189–203 (1999).
    https://doi.org/10.1007/s006070050002.
"""
struct DelCorsoManzini <: AbstractDecider end

push!(ALGORITHMS[:Recognition], DelCorsoManzini)

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
tracking the latest positions at which the remaining indices can be placed [DM99].

The incorporation of perimeter search to this approach entails precomputing a "perimeter" of
``d``-permutations of row indices of ``A``, where ``d`` is a positive integer passed as a
parameter to the decider. Each permutation represents a way to select the last ``d`` entries
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
- `depth::D<:Union{Nothing,Integer}`: the perimeter search depth. If this field is not set
    (and thus automatically initialized to `nothing`), a default depth is computed by
    [`dcm_ps_optimal_depth`](@ref) as a function of the input matrix every time the decider
    is passed to [`has_bandwidth_k_ordering`](@ref) as a function of the input matrix.
    Otherwise, it must be manually set to a positive integer.

# Constructors
- `DelCorsoManziniWithPS()`: constructs a new `DelCorsoManziniWithPS` instance with the
    default perimeter search depth initialized to `nothing`.
- `DelCorsoManziniWithPS(depth::Integer)`: constructs a new `DelCorsoManziniWithPS` instance
    with the specified perimeter search depth. `depth` must be a positive integer.

# Supertype Hierarchy
`DelCorsoManziniWithPS` <: [`AbstractDecider`](@ref) <: [`AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix ``A``, perimeter search depth ``d``, and threshold bandwidth
``k``, the Del Corso–Manzini algorithm with perimeter search runs in ``O(n! ⋅ max(nᵈ, nk))``
time:
- We perform a depth-first search of ``O(n!)`` partial orderings.
- Checking plausibility of each partial ordering takes ``O(nk)`` time, and checking
    compatibility with all size-``d`` LPOs takes ``O(nᵈ)`` time. Thus, the overall time
    complexity for each value of ``k`` is ``O(n! ⋅ (nᵈ + nk))``.
- Therefore, the overall time complexity is ``O(n! ⋅ max(nᵈ, nk))``.

Of course, this is but an upper bound on the time complexity of Del Corso–Manzini with
perimeter search, achieved only in the most pathological of cases. In practice, efficient
pruning techniques and compatibility checks result in approximately exponential growth in
time complexity with respect to ``n``.

Based on experimental results, the algorithm is feasible for ``n×n`` matrices up to
``n ≈ 100`` or so.

# Examples
Here, Del Corso–Manzini with perimeter search ascertains that A random ``30×30`` matrix has
a minimum bandwidth greater than ``9``. The depth parameter is not explicitly set; instead,
some near-optimal value is automatically computed upon the first
[`has_bandwidth_k_ordering`](@ref) function call.
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(5847);

julia> (n, p) = (30, 0.05);

julia> A = sprand(n, n, p);

julia> A = A + A' # Ensure structural symmetry;

julia> k = 6;

julia> has_bandwidth_k_ordering(A, k, Recognition.DelCorsoManziniWithPS())
Results of Bandwidth Recognition Algorithm
 * Algorithm: Del Corso–Manzini with perimeter search
 * Bandwidth Threshold k: 6
 * Has Bandwidth ≤ k Ordering: false
 * Original Bandwidth: 27
 * Matrix Size: 30×30
```

Now, Del Corso–Manzini with perimeter search recognizes that a random ``35×35`` matrix has a
minimum bandwidth at most ``8``. In this case, we explitily set the depth parameter to ``4``
beforehand instead of relying on [`Recognition.dcm_ps_optimal_depth`](@ref).
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(23552);

julia> (n, p, depth) = (35, 0.02, 4);

julia> A = sprand(n, n, p);

julia> A = A + A' # Ensure structural symmetry;

julia> k = 8;

julia> has_bandwidth_k_ordering(A, k, Recognition.DelCorsoManziniWithPS(depth))
Results of Bandwidth Recognition Algorithm
 * Algorithm: Del Corso–Manzini with perimeter search
 * Bandwidth Threshold k: 8
 * Has Bandwidth ≤ k Ordering: true
 * Original Bandwidth: 32
 * Matrix Size: 35×35
```

# Notes
For readers of the original paper, what we call the Del Corso–Manzini recognition algorithm
with perimeter search here is essentially a wrapper around the underlying `AddNode1` and
`Prune` subroutines in what [DM99, p. 193] term the "MB-PS algorithm" for bandwidth
minimization (not mere recognition). MB-PS (which we also implement in
[`MatrixBandwidth.Minimization.Exact.DelCorsoManziniWithPS`](@ref)) calls this recognition
procedure with incrementing values of ``k`` until a bandwidth-``k`` ordering is found, with
``k`` initialized to some lower bound on the minimum bandwidth of ``A`` up to symmetric
permutation.

[DM99, p. 191] also describe an "MB-ID algorithm" for bandwidth minimization, which we
implement in [`MatrixBandwidth.Minimization.Exact.DelCorsoManzini`](@ref). Similarly, the
underlying recognition subroutine for MB-ID is implemented in [`DelCorsoManzini`](@ref).

# References
- [DM99](@cite): G. M. Del Corso and G. Manzini. *Finding Exact Solutions to the Bandwidth
    Minimization Problem*. Computing **62**, 189–203 (1999).
    https://doi.org/10.1007/s006070050002.
"""
struct DelCorsoManziniWithPS{D<:Union{Nothing,Integer}} <: AbstractDecider
    depth::D

    #= We cannot compute a (hopefully) near-optimal perimeter search depth upon
    instantiation of the decider, as it depends on the input matrix as well. Hence, we use
    `nothing` as a sentinel to indicate to `_bool_bandwidth_k_ordering` that a default depth
    still needs to be computed upon the function call. =#
    DelCorsoManziniWithPS() = new{Nothing}(nothing)

    function DelCorsoManziniWithPS(depth::Integer)
        if depth <= 0
            throw(ArgumentError("Perimeter search depth must be positive, got $depth"))
        end

        return new{Integer}(depth)
    end
end

push!(ALGORITHMS[:Recognition], DelCorsoManziniWithPS)

Base.summary(::DelCorsoManziniWithPS) = "Del Corso–Manzini with perimeter search"

_requires_symmetry(::DelCorsoManziniWithPS) = true

"""
    dcm_ps_optimal_depth(A) -> Int

Compute a (hopefully) near-optimal Del Corso–Manzini perimeter search depth for `A`.

Taking experimental results from [DM99, pp. 197–99] into account, this function tries to
approximate the optimal depth parameter as a function of both matrix size and density. This
depth parameter determines how large of a "perimeter" of last-placed indices is precomputed
in the Del Corso–Manzini algorithm with perimeter search.

# Arguments
- `A::AbstractMatrix{Bool}`: the (structurally symmetric and square) input matrix whose
    bandwidth is investigated.

# Returns
- `::Int`: a (hopefully) near-optimal perimeter search depth for the Del Corso–Manzini
    algorithm with perimeter search on `A`.

# Notes
See Tables 4, 5, and 6 from the original paper for more details on experimental results
regarding the optimal perimeter search depth [DM99, pp. 197–99].

See also [`DelCorsoManziniWithPS`](@ref) and
[`MatrixBandwidth.Minimization.Exact.DelCorsoManziniWithPS`](@ref)) for our implementation
of the relevant bandwidth recognition and bandwidth minimization algorithms, respectively.

# References
- [DM99](@cite): G. M. Del Corso and G. Manzini. *Finding Exact Solutions to the Bandwidth
    Minimization Problem*. Computing **62**, 189–203 (1999).
    https://doi.org/10.1007/s006070050002.
"""
function dcm_ps_optimal_depth(A::AbstractMatrix{Bool})
    n = size(A, 1)

    #= `A` is symmetric and diagonal entries are irrelevant (and indeed all false in this
    case), so we need only count entries above the main diagonal. =#
    num_nonzero = 2count(A[i, j] for i in 1:(n - 1) for j in (i + 1):n)

    #= Compute a lower bound on the bandwidth up to symmetric permutation and hence an upper
    bound on the band density given a minimal bandwidth ordering. =#
    k_lb = bandwidth_lower_bound(A)
    band_density_ub = num_nonzero / (2n * k_lb - k_lb * (k_lb + 1))

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

function _bool_bandwidth_k_ordering(A::AbstractMatrix{Bool}, k::Integer, ::DelCorsoManzini)
    n = size(A, 1)

    ordering_buf = Vector{Int}(undef, n)
    adj_lists = map(node -> findall(view(A, :, node)), 1:n)

    unselected = Set(1:n)
    adj_list = Set{Int}()
    #= Sentinel value for a nonempty perimeter just so the common logic between DCM and
    DCM-PS still runs. =#
    perimeter = [(Int[], Int[])]
    num_placed = 0
    ps_depth = 0

    return _dcm_add_node!(
        ordering_buf, A, k, adj_lists, unselected, adj_list, perimeter, num_placed, ps_depth
    )
end

function _bool_bandwidth_k_ordering(
    A::AbstractMatrix{Bool}, k::Integer, ::DelCorsoManziniWithPS{Nothing}
)
    #= We cannot compute a (hopefully) near-optimal perimeter search depth upon
    instantiation of the decider, as it depends on the input matrix as well. =#
    decider_with_depth = DelCorsoManziniWithPS(dcm_ps_optimal_depth(A))
    return _bool_bandwidth_k_ordering(A, k, decider_with_depth)
end

function _bool_bandwidth_k_ordering(
    A::AbstractMatrix{Bool}, k::Integer, decider::DelCorsoManziniWithPS{Integer}
)
    n = size(A, 1)
    ps_depth = decider.depth

    if ps_depth > n
        throw(ArgumentError("Perimeter search depth $ps_depth exceeds matrix order $n"))
    end

    ordering_buf = Vector{Int}(undef, n)
    adj_lists = map(node -> findall(view(A, :, node)), 1:n)

    unselected = Set(1:n)
    adj_list = Set{Int}()
    num_placed = 0
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
    k::Integer,
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
    indices without violating the bandwidth-`k` constraint.

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
        far_nodes = view(ordering_buf, 1:(num_placed - k))

        #= Ensure that no placed node at least `k` positions before `candidate` is adjacent
        to `candidate` (otherwise, the bandwidth would be at least `k + 1`). =#
        if !any(view(A, far_nodes, candidate))
            unselected_new = setdiff(unselected, [candidate])

            if !_dcm_order_is_reversed(
                ordering_buf, unselected_new, num_placed, ps_depth, candidate
            )
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
        end

        res = iterate(unselected, state)
    end

    return ordering
end

#= `i₁, i₂, … iₙ` induces the same bandwidth as `iₙ, iₙ₋₁, … i₁`, so without loss of
generality, we restrict our search to partial orderings such that `i₁` is less than or equal
to the maximum possible `iₙ` (with equality checked just in case `n = 1`). =#
function _dcm_order_is_reversed(
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
    placed at an intermediate position that we can ignore. (Also, we avoid reducing over an
    empty collection with `maximum` and throwing an error…) =#
    if isempty(unselected)
        max_last_label = candidate
    else
        max_last_label = maximum(unselected)
    end

    return max_last_label <= first_label
end

function _dcm_lpo_time_stamps(lpo::Vector{Int}, A::AbstractMatrix{Bool}, k::Integer)
    n = size(A, 1)

    time_stamps = zeros(Int, n)
    d = length(lpo)
    #= The nodes already in the LPO occupy fixed positions (in the last `d`), so we have
    trivial tight lower bounds on the earliest positions at which they can be placed (and
    indeed are precisely placed in this ordering). =#
    foreach(((i, node),) -> time_stamps[node] = n - d + i, enumerate(lpo))

    queue = Queue{Int}()
    foreach(node -> enqueue!(queue, node), lpo)

    #= The nodes processed here are those which are not fixed in the last `d` positions, so
    we compute loose lower bounds on the earliest positions at which they can be placed. =#
    while !isempty(queue)
        node = dequeue!(queue)
        unvisited = filter!(
            neighbor -> time_stamps[neighbor] == 0, findall(view(A, :, node))
        )
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

        #= `i₁, i₂, … iₙ` induces the same bandwidth as `iₙ, iₙ₋₁, … i₁`, so without loss of
        generality, we restrict our search to orderings such that `i₁ ≤ iₙ` (with equality
        checked just in case `n = 1`). =#
        if num_placed == 0
            filter!(((lpo, _),) -> candidate <= lpo[end], perimeter_new)
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

#= With the usual bandwidth-`k` constraint for already placed nodes checked inline in
`_dcm_add_node!`, this implements the expiration-time pruning (for both DCM and DCM-PS). =#
function _dcm_is_compatible(
    ordering_buf::Vector{Int},
    A::AbstractMatrix{Bool},
    adj_list::Set{Int},
    k::Integer,
    num_placed::Int,
)
    l = length(adj_list)

    if l > k
        #= If more than `k` nodes are adjacent to any node (not just the most recent one)
        already placed, then it is impossible to place them all from `num_placed + 1` to
        `num_placed + k`, thus forcing a violation of the bandwidth-`k` constraint. =#
        #= If more than `k` unplaced nodes are adjacent to our candidate, we cannot place
        them all after said candidate without violating the bandwidth-`k` constraint. =#
        compat = false
    else
        placed_nodes = view(ordering_buf, 1:num_placed)
        constraints = (1:l) .+ num_placed

        latest_positions = Iterators.map(
            neighbor -> findfirst(node -> A[node, neighbor], placed_nodes), adj_list
        )
        #= We must add `k` to the latest allowed positions to account for the bandwidth
        constraint. To handle cases where some `neighbor` is only adjacent to the candidate
        and not to any already placed `node` (and thus the `findfirst` call returns
        `nothing`), we use the sentinel value `typemax(Int)`. =#
        latest_positions = Iterators.map(
            pos -> isnothing(pos) ? typemax(Int) : pos + k, latest_positions
        )
        latest_positions = sort!(collect(latest_positions))

        compat = all(latest_positions .>= constraints)
    end

    return compat
end
