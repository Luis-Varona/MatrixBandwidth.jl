# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    DelCorsoManzini <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

The *Del Corso–Manzini minimization algorithm* is an exact method for minimizing the
bandwidth of a structurally symmetric matrix ``A``. For a fixed ``k ∈ ℕ``, the algorithm
performs a depth-first search of all partial orderings of the rows and columns of ``A``,
adding indices one at a time. Partial orderings are pruned not only by ensuring that
adjacent pairs of currently placed indices are within ``k`` of each other but also by
tracking the latest positions at which the remaining indices can be placed. This search is
repeated with incrementing values of ``k`` until a bandwidth-``k`` ordering is found
[DCM99](@cite), with ``k`` initialized to some lower bound on the minimum bandwidth of ``A``
up to symmetric permutation.

Specifically, this implementation of the Del Corso–Manzini algorithm uses the
``min(α(A), γ(A))`` lower bound from [CSG05; pp. 359--60](@cite) as the initial value of
``k``. (Further implementation details can be found in the source code for
[`bandwidth_lower_bound`](@ref).) This improves upon the original algorithm, which used the
maximum number of nonzero off-diagonal entries in a single row as a lower bound on the
minimum bandwidth of ``A`` up to symmetric permutation [DCM99; p. 192--93](@cite).

As noted above, the Del Corso–Manzini algorithm requires structurally symmetric input (that
is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for ``1 ≤ i, j ≤ n``).

# Performance
Given an ``n×n`` input matrix ``A``, the Del Corso–Manzini algorithm runs in
``O(n³ ⋅ n!)`` time:
- For each underlying "bandwidth ≤ ``k``" check, we perform a depth-first search of
    ``O(n!)`` partial orderings.
- Checking plausibility of each partial ordering takes ``O(nk)`` time, resulting in
    ``O(nk ⋅ n!)`` steps for each value of ``k``.
- The difference between the maximum possible bandwidth (``n - 1``) and our initial lower
    bound grows linearly in ``n``, so we run the underlying ``O(nk ⋅ n!)`` recognition
    algorithm ``O(n)`` times.
- Finally, ``∑ₖ₌₀ⁿ⁻¹ nk = O(n³)``, so the overall time complexity is ``O(n³ ⋅ n!)``.

Of course, this is but an upper bound on the time complexity of Del Corso–Manzini, achieved
only in the most pathological of cases. In practice, efficient pruning techniques and
compatibility checks—along with [CSG05; pp. 359--60](@cite)'s relatively tight initial lower
bound on the minimum bandwidth—result in approximately exponential growth in time complexity
with respect to ``n``.

Based on experimental results, the algorithm is feasible for ``n×n`` matrices up to
``n ≈ 100`` or so.

# Examples
[TODO: Write here]

# Notes
For readers of the original paper, what we call the Del Corso–Manzini minimization algorithm
here is designated the "MB-ID algorithm" in [DCM99; p. 191](@cite). The so-called "MB-PS
algorithm," on the other hand, we implement in [`Minimization.DelCorsoManziniWithPS`](@ref).
"""
struct DelCorsoManzini <: ExactSolver end

Base.summary(::DelCorsoManzini) = "Del Corso–Manzini"

_requires_symmetry(::DelCorsoManzini) = true

"""
    DelCorsoManziniWithPS{D} <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

The *Del Corso–Manzini minimization algorithm with perimeter search* is an exact method for
minimizing the bandwidth of a structurally symmetric matrix ``A``. The base Del
Corso–Manzini algorithm performs a depth-first search of all partial orderings of the rows
and columns of ``A`` for some fixed ``k ∈ ℕ``, adding indices one at a time. Partial
orderings are pruned not only by ensuring that adjacent pairs of currently placed indices
are within ``k`` of each other but also by tracking the latest positions at which the
remaining indices can be placed. This search is repeated with incrementing values of ``k``
until a bandwidth-``k`` ordering is found [DCM99](@cite), with ``k`` initialized to some
lower bound on the minimum bandwidth of ``A`` up to symmetric permutation.

The incorporation of perimeter search to this approach entails precomputing a "perimeter" of
``d``-permutations of row indices of ``A``, where ``d`` is a positive integer passed as a
parameter to the solver. Each permutation represents a way to select the last ``d`` entries
of the ordering, and as the construction of the partial ordering progresses, potential
endings are pruned to exclude those incompatible with already placed indices. In addition to
pruning a potential ending if it contains indices already placed, compatibility is also
checked via precomputed time stamps indicating, for each potential ending, a loose lower
bound on the earliest position at which any given index can be placed should said ending be
selected.

Like our implementation of the base Del Corso–Manzini algorithm (see
[`Minimization.DelCorsoManzini`](@ref)), this implementation uses the
``min(α(A), γ(A))`` lower bound from [CSG05; pp. 359--60](@cite) as the initial value of
``k``. (Further implementation details can be found in the source code for
[`bandwidth_lower_bound`](@ref).) This improves upon the original algorithm, which used the
maximum number of nonzero off-diagonal entries in a single row as a lower bound on the
minimum bandwidth of ``A`` up to symmetric permutation [DCM99; p. 194](@cite).

As noted above, the Del Corso–Manzini algorithm with perimeter search requires structurally
symmetric input (that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is
nonzero for ``1 ≤ i, j ≤ n``).

# Fields
[TODO: Write here]

# Constructors
[TODO: Write here]

# Performance
Given an ``n×n`` input matrix ``A`` and perimeter search depth ``d``, the Del Corso–Manzini
algorithm with perimeter search runs in ``O(nᴰ⁺¹ ⋅ n!)`` time, where ``Dᴰ = max(d, 2)``:
- For each underlying "bandwidth ≤ ``k``" check, we perform a depth-first search of
    ``O(n!)`` partial orderings.
- Checking plausibility of each partial ordering takes ``O(nk)`` time, and checking
    compatibility with all size-``d`` LPOs takes ``O(nᵈ)`` time. Thus, the overall time
    complexity for each value of ``k`` is ``O((nᵈ + nk) ⋅ n!)``.
- The difference between the maximum possible bandwidth (``n - 1``) and our initial lower
    bound grows linearly in ``n``, so we run the underlying ``O((nᵈ + nk) ⋅ n!)``
    recognition algorithm ``O(n)`` times.
- Finally, ``∑ₖ₌₀ⁿ⁻¹ (nᵈ + nk) = O(nᵈ⁺¹ + n³)``, so the overall time complexity
    is ``O(nᴰ⁺¹ ⋅ n!)``, where ``D = max(d, 2)``.

Of course, this is but an upper bound on the time complexity of Del Corso–Manzini with
perimeter search, achieved only in the most pathological of cases. In practice, efficient
pruning techniques and compatibility checks—along with [CSG05; pp. 359--60](@cite)'s
relatively tight initial lower bound on the minimum bandwidth—result in approximately
exponential growth in time complexity with respect to ``n``.

Based on experimental results, the algorithm is feasible for ``n×n`` matrices up to
``n ≈ 100`` or so.

# Examples
[TODO: Write here]

# Notes
For readers of the original paper, what we call the Del Corso–Manzini minimization algorithm
with perimeter search here is designated the "MB-PS algorithm" in [DCM99; p. 193](@cite).
The so-called "MB-ID algorithm," on the other hand, we implement in
[`Minimization.DelCorsoManzini`](@ref).
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

# TODO: Add inline comments to all the code in the rest of the file below

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
        perimeter = map(lpo -> (lpo, Recognition._dcm_lpo_time_stamps(lpo, A, k)), lpos)
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
