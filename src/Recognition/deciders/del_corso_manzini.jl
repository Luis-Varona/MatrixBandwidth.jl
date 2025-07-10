# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    DelCorsoManzini <: AbstractDecider <: AbstractAlgorithm

TODO: Write here
"""
struct DelCorsoManzini <: AbstractDecider end

Base.summary(::DelCorsoManzini) = "Del Corso–Manzini"

_requires_symmetry(::DelCorsoManzini) = true

"""
    DelCorsoManziniWithPS{D} <: AbstractDecider <: AbstractAlgorithm

TODO: Write here
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

# TODO: Comment here
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
        (lpo, _lpo_time_stamps(lpo, A, k)) for node_subset in combinations(1:n, ps_depth)
        for lpo in permutations(node_subset)
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

            #= If `num_placed` is zero and we check `ordering_buf[1] < max_last_label`, we
            risk erroneous behavior if `ordering_buf[1]` is initialized to a value greater
            than or equal to `max_last_label` (e.g., when using `UndefInitializer()`). =#
            if num_placed == 0
                first_label = candidate
            else
                first_label = ordering_buf[1]
            end

            # `max_label` is maximum possible label of the last node that may be placed
            if isempty(unselected_new)
                max_last_label = candidate
            else
                max_last_label = maximum(unselected_new)
            end

            #= The ordering `i₁, i₂, … iₙ` induces the same bandwidth as `iₙ, iₙ₋₁, … i₁`, so
            without loss of generality, we restrict our search to partial orderings such
            that `i₁` is less than the maximum possible `iₙ`. =#
            if first_label < max_last_label
                perimeter_new = _dcm_pruned_perimeter(
                    perimeter, candidate, num_placed, ps_depth
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

function _dcm_pruned_perimeter(
    perimeter::Vector{Tuple{Vector{Int},Vector{Int}}},
    candidate::Int,
    num_placed::Int,
    search_depth::Int,
)
    if search_depth == 0
        perimeter_new = perimeter
    else
        perimeter_new = filter(
            ((_, time_stamp),) -> time_stamp[candidate] <= num_placed + 1, perimeter
        )

        if num_placed == 0
            filter!(((lpo, _),) -> lpo[search_depth] > candidate, perimeter_new)
        end
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

function _lpo_time_stamps(lpo::Vector{Int}, A::AbstractMatrix{Bool}, k::Int)
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
