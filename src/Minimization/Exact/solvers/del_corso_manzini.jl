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
[DM99], with ``k`` initialized to some lower bound on the minimum bandwidth of ``A`` up to
symmetric permutation.

Specifically, this implementation of the Del Corso–Manzini algorithm uses the
``min(α(A), γ(A))`` lower bound from [CS05, pp. 359--60] as the initial value of ``k``.
(Further implementation details can be found in the source code for
[`bandwidth_lower_bound`](@ref).) This improves upon the original algorithm, which used the
maximum number of nonzero off-diagonal entries in a single row as a lower bound on the
minimum bandwidth of ``A`` up to symmetric permutation [DM99, p. 192--93].

As noted above, the Del Corso–Manzini algorithm requires structurally symmetric input (that
is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for ``1 ≤ i, j ≤ n``).

# Supertype Hierarchy
`DelCorsoManzini` <: [`ExactSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix ``A``, the Del Corso–Manzini algorithm runs in
``O(n! ⋅ n³)`` time:
- For each underlying "bandwidth ≤ ``k``" check, we perform a depth-first search of
    ``O(n!)`` partial orderings.
- Checking plausibility of each partial ordering takes ``O(nk)`` time, resulting in
    ``O(n! ⋅ nk)`` steps for each value of ``k``.
- The difference between the maximum possible bandwidth (``n - 1``) and our initial lower
    bound grows linearly in ``n``, so we run the underlying ``O(n! ⋅ nk)`` recognition
    algorithm ``O(n)`` times.
- Finally, ``∑ₖ₌₀ⁿ⁻¹ nk = O(n³)``, so the overall time complexity is ``O(n! ⋅ n³)``.

Of course, this is but an upper bound on the time complexity of Del Corso–Manzini, achieved
only in the most pathological of cases. In practice, efficient pruning techniques and
compatibility checks—along with [CS05, pp. 359--60]'s relatively tight initial lower bound
on the minimum bandwidth—result in approximately exponential growth in time complexity with
respect to ``n``.

Based on experimental results, the algorithm is feasible for ``n×n`` matrices up to
``n ≈ 100`` or so.

# Examples
We verify the optimality of the ordering found by Del Corso–Manzini for a random ``9×9``
matrix via a brute-force search over all possible permutations up to reversal:
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(0117);

julia> (n, p) = (9, 0.5);

julia> A = sprand(n, n, p);

julia> A = A + A' # Ensure structural symmetry;

julia> res_bf = minimize_bandwidth(A, Minimization.BruteForceSearch())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Brute-force search
 * Approach: exact
 * Minimum Bandwidth: 5
 * Original Bandwidth: 8
 * Matrix Size: 9×9

julia> res_dcm = minimize_bandwidth(A, Minimization.DelCorsoManzini())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Del Corso–Manzini
 * Approach: exact
 * Minimum Bandwidth: 5
 * Original Bandwidth: 8
 * Matrix Size: 9×9
```

We now generate (and shuffle) a random ``40×40`` matrix with minimum bandwidth ``10`` using
[`MatrixBandwidth.random_banded_matrix`](@ref). Del Corso–Manzini then finds a
bandwidth-``10`` ordering, which is (we claim) optimal up to symmetric permutation. (In some
cases, `random_banded_matrix(n, k)` *does* generate matrices with minimum bandwidth `< k`.
Nevertheless, this example demonstrates that Del Corso–Manzini at the very least finds a
good ordering, even though exact optimality—which *is* guaranteed by the original paper
[DM99]—is not explicitly verified.)
```jldoctest
julia> using Random

julia> Random.seed!(0201);

julia> (n, k) = (40, 10);

julia> A = random_banded_matrix(n, k);

julia> perm = randperm(n);

julia> A_shuffled = A[perm, perm];

julia> bandwidth(A)
10

julia> bandwidth(A_shuffled) # Much larger after shuffling
36

julia> minimize_bandwidth(A_shuffled, Minimization.DelCorsoManzini())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Del Corso–Manzini
 * Approach: exact
 * Minimum Bandwidth: 10
 * Original Bandwidth: 36
 * Matrix Size: 40×40
```

# Notes
For readers of the original paper, what we call the Del Corso–Manzini minimization algorithm
here is designated the "MB-ID algorithm" in [DM99, p. 191]. The so-called "MB-PS algorithm,"
on the other hand, we implement in [`DelCorsoManziniWithPS`](@ref).

# References

- [CS05](@cite): A. Caprara and J.-J. Salazar-González. *Laying Out Sparse Graphs with
    Provably Minimum Bandwidth*. INFORMS Journal on Computing **17**, 356–73 (2005).
    https://doi.org/10.1287/ijoc.1040.0083.
- [DM99](@cite): G. M. Del Corso and G. Manzini. *Finding Exact Solutions to the Bandwidth
    Minimization Problem*. Computing **62**, 189–203 (1999).
    https://doi.org/10.1007/s006070050002.
"""
struct DelCorsoManzini <: ExactSolver end

push!(ALGORITHMS[:Minimization][:Exact], DelCorsoManzini)

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
until a bandwidth-``k`` ordering is found [DM99], with ``k`` initialized to some lower bound
on the minimum bandwidth of ``A`` up to symmetric permutation.

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
[`DelCorsoManzini`](@ref)), this implementation uses the ``min(α(A), γ(A))`` lower bound
from [CS05, pp. 359--60] as the initial value of ``k``. (Further implementation details can
be found in the source code for [`bandwidth_lower_bound`](@ref).) This improves upon the
original algorithm, which used the maximum number of nonzero off-diagonal entries in a
single row as a lower bound on the minimum bandwidth of ``A`` up to symmetric permutation
[DM99, p. 194].

As noted above, the Del Corso–Manzini algorithm with perimeter search requires structurally
symmetric input (that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is
nonzero for ``1 ≤ i, j ≤ n``).

# Fields
- `depth::D<:Union{Nothing,Int}`: the perimeter search depth. If this field is not set (and)
    thus automatically initialized to `nothing`), a default depth is automatically computed
    by [`Recognition.dcm_ps_optimal_depth`](@ref) as a function of the input matrix every
    time the solver is passed to [`MatrixBandwidth.Minimization.minimize_bandwidth`](@ref).
    Otherwise, it must be manually set to a positive integer.

# Constructors
- `DelCorsoManziniWithPS()`: constructs a new `DelCorsoManziniWithPS` instance with the
    default perimeter search depth initialized to `nothing`.
- `DelCorsoManziniWithPS(depth::Integer)`: constructs a new `DelCorsoManziniWithPS` instance
    with the specified perimeter search depth. `depth` must be a positive integer.

# Supertype Hierarchy
`DelCorsoManziniWithPS` <: [`ExactSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix ``A`` and perimeter search depth ``d``, the Del Corso–Manzini
algorithm with perimeter search runs in ``O(n! ⋅ nᴰ⁺¹)`` time, where ``Dᴰ = max(d, 2)``:
- For each underlying "bandwidth ≤ ``k``" check, we perform a depth-first search of
    ``O(n!)`` partial orderings.
- Checking plausibility of each partial ordering takes ``O(nk)`` time, and checking
    compatibility with all size-``d`` LPOs takes ``O(nᵈ)`` time. Thus, the overall time
    complexity for each value of ``k`` is ``O(n! ⋅ (nᵈ + nk))``.
- The difference between the maximum possible bandwidth (``n - 1``) and our initial lower
    bound grows linearly in ``n``, so we run the underlying ``O(n! ⋅ (nᵈ + nk))``
    recognition algorithm ``O(n)`` times.
- Finally, ``∑ₖ₌₀ⁿ⁻¹ (nᵈ + nk) = O(nᵈ⁺¹ + n³)``, so the overall time complexity
    is ``O(n! ⋅ nᴰ⁺¹)``, where ``D = max(d, 2)``.

Of course, this is but an upper bound on the time complexity of Del Corso–Manzini with
perimeter search, achieved only in the most pathological of cases. In practice, efficient
pruning techniques and compatibility checks—along with [CS05, pp. 359--60]'s relatively
tight initial lower bound on the minimum bandwidth—result in approximately exponential
growth in time complexity with respect to ``n``.

Based on experimental results, the algorithm is feasible for ``n×n`` matrices up to
``n ≈ 100`` or so.

# Examples
We verify the optimality of the ordering found by Del Corso–Manzini with perimeter search
for a random ``9×9`` matrix via a brute-force search over all possible permutations up to
reversal. The depth parameter is not explicitly set; instead, some near-optimal value is
automatically computed upon the first
[`MatrixBandwidth.Minimization.minimize_bandwidth`](@ref) function call.
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(548836);

julia> (n, p) = (9, 0.2);

julia> A = sprand(n, n, p);

julia> A = A + A' # Ensure structural symmetry;

julia> res_bf = minimize_bandwidth(A, Minimization.BruteForceSearch())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Brute-force search
 * Approach: exact
 * Minimum Bandwidth: 3
 * Original Bandwidth: 8
 * Matrix Size: 9×9

julia> res_dcm = minimize_bandwidth(A, Minimization.DelCorsoManziniWithPS())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Del Corso–Manzini with perimeter search
 * Approach: exact
 * Minimum Bandwidth: 3
 * Original Bandwidth: 8
 * Matrix Size: 9×9
```

We now generate (and shuffle) a random ``30×30`` matrix with minimum bandwidth ``8`` using
[`MatrixBandwidth.random_banded_matrix`](@ref). Del Corso–Manzini with perimeter search then
finds a bandwidth-``8`` ordering, which is (we claim) optimal up to symmetric permutation.
(In some cases, `random_banded_matrix(n, k)` *does* generate matrices with minimum bandwidth
`< k`. Nevertheless, this example demonstrates that Del Corso–Manzini at the very least
finds a good ordering, even though exact optimality—which *is* guaranteed by the original
paper [DM99]—is not explicitly verified.) In this case, we set the depth parameter to ``4``
beforehand instead of relying on [`Recognition.dcm_ps_optimal_depth`](@ref).

```jldoctest
julia> using Random

julia> Random.seed!(78779);

julia> (n, k, depth) = (30, 8, 4);

julia> A = random_banded_matrix(n, k);

julia> perm = randperm(n);

julia> A_shuffled = A[perm, perm];

julia> bandwidth(A)
8

julia> bandwidth(A_shuffled) # Much larger after shuffling
25

julia> minimize_bandwidth(A_shuffled, Minimization.DelCorsoManziniWithPS(depth))
Results of Bandwidth Minimization Algorithm
 * Algorithm: Del Corso–Manzini with perimeter search
 * Approach: exact
 * Minimum Bandwidth: 8
 * Original Bandwidth: 25
 * Matrix Size: 30×30
```

# Notes
For readers of the original paper, what we call the Del Corso–Manzini minimization algorithm
with perimeter search here is designated the "MB-PS algorithm" in [DM99. p. 193]. The
so-called "MB-ID algorithm," on the other hand, we implement in [`DelCorsoManzini`](@ref).

# References

- [CS05](@cite): A. Caprara and J.-J. Salazar-González. *Laying Out Sparse Graphs with
    Provably Minimum Bandwidth*. INFORMS Journal on Computing **17**, 356–73 (2005).
    https://doi.org/10.1287/ijoc.1040.0083.
- [DM99](@cite): G. M. Del Corso and G. Manzini. *Finding Exact Solutions to the Bandwidth
    Minimization Problem*. Computing **62**, 189–203 (1999).
    https://doi.org/10.1007/s006070050002.
"""
struct DelCorsoManziniWithPS{D<:Union{Nothing,Integer}} <: ExactSolver
    depth::D

    #= We cannot compute a (hopefully) near-optimal perimeter search depth upon
    instantiation of the solver, as it depends on the input matrix as well. Hence, we use
    `nothing` as a sentinel to indicate to `_bool_minimal_band_ordering` that a default
    depth still needs to be computed upon the function call. =#
    DelCorsoManziniWithPS() = new{Nothing}(nothing)

    function DelCorsoManziniWithPS(depth::Integer)
        if depth <= 0
            throw(ArgumentError("Perimeter search depth must be positive, got $depth"))
        end

        return new{Integer}(depth)
    end
end

push!(ALGORITHMS[:Minimization][:Exact], DelCorsoManziniWithPS)

Base.summary(::DelCorsoManziniWithPS) = "Del Corso–Manzini with perimeter search"

_requires_symmetry(::DelCorsoManziniWithPS) = true

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::DelCorsoManzini)
    n = size(A, 1)

    ordering_buf = Vector{Int}(undef, n)
    k = bandwidth_lower_bound(A)
    adj_lists = map(node -> findall(view(A, :, node)), 1:n)

    unselected = Set(1:n)
    adj_list = Set{Int}()
    #= Sentinel value for a nonempty perimeter just so the common logic between DCM and
    DCM-PS still runs. =#
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
    #= We cannot compute a (hopefully) near-optimal perimeter search depth upon
    instantiation of the solver, as it depends on the input matrix as well. =#
    solver_with_depth = DelCorsoManziniWithPS(Recognition.dcm_ps_optimal_depth(A))
    return _bool_minimal_band_ordering(A, solver_with_depth)
end

function _bool_minimal_band_ordering(
    A::AbstractMatrix{Bool}, solver::DelCorsoManziniWithPS{Integer}
)
    n = size(A, 1)
    ps_depth = solver.depth

    if ps_depth > n
        throw(ArgumentError("Perimeter search depth $ps_depth exceeds matrix order $n"))
    end

    ordering_buf = Vector{Int}(undef, n)
    k = bandwidth_lower_bound(A)
    adj_lists = map(node -> findall(view(A, :, node)), 1:n)

    unselected = Set(1:n)
    adj_list = Set{Int}()
    num_placed = 0
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
