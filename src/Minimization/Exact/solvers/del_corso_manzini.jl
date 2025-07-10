# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    DelCorsoManzini <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

TODO: Write here
"""
struct DelCorsoManzini <: ExactSolver end

Base.summary(::DelCorsoManzini) = "Del Corso–Manzini"

_requires_symmetry(::DelCorsoManzini) = true

"""
    DelCorsoManziniWithPS{D} <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

TODO: Write here
"""
struct DelCorsoManziniWithPS{D<:Union{Nothing,Int}} <: ExactSolver
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

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::DelCorsoManzini)
    n = size(A, 1)

    ordering_buf = Vector{Int}(undef, n)
    k = bandwidth_lower_bound(A)
    adj_lists = map(node -> findall(A[:, node]), 1:n)

    unselected = Set(1:n)
    adj_list = Set{Int}()
    perimeter = [(Int[], Int[])]
    num_placed = 0
    ps_depth = 0

    ordering = nothing

    while isnothing(ordering)
        ordering = Recognition._dcm_add_node!(
            ordering_buf,
            A,
            k,
            adj_lists,
            unselected,
            adj_list,
            perimeter,
            num_placed,
            ps_depth,
        )
        k += 1
    end

    return ordering
end

function _bool_minimal_band_ordering(
    A::AbstractMatrix{Bool}, solver::DelCorsoManziniWithPS{Nothing}
)
    decider_with_depth = DelCorsoManziniWithPS(Recognition._dcm_ps_optimal_depth(A))
    return _bool_minimal_band_ordering(A, decider_with_depth)
end

function _bool_minimal_band_ordering(
    A::AbstractMatrix{Bool}, solver::DelCorsoManziniWithPS{Int}
)
    n = size(A, 1)

    ordering_buf = Vector{Int}(undef, n)
    k = bandwidth_lower_bound(A)
    adj_lists = map(node -> findall(A[:, node]), 1:n)

    unselected = Set(1:n)
    adj_list = Set{Int}()
    num_placed = 0
    ps_depth = solver.depth
    lpos = Iterators.flatmap(permutations, combinations(1:n, ps_depth))

    ordering = nothing

    while isnothing(ordering)
        perimeter = map(lpo -> (lpo, Recognition._lpo_time_stamps(lpo, A, k)), lpos)
        ordering = Recognition._dcm_add_node!(
            ordering_buf,
            A,
            k,
            adj_lists,
            unselected,
            adj_list,
            perimeter,
            num_placed,
            ps_depth,
        )
        k += 1
    end

    return ordering
end
