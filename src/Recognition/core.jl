# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    has_bandwidth_k_ordering(A, k, decider=CapraraSalazarGonzalez()) -> RecognitionResult

Determine whether `A` has bandwidth at most `k` using the algorithm defined by `decider`.

The *bandwidth* of an ``n×n`` matrix ``A`` is the minimum non-negative integer
``k ∈ \\{0, 1, …, n - 1\\}`` such that ``A[i, j] = 0`` whenever ``|i - j| > k``.

Given some fixed non-negative integer `k`, this function determines (with 100% certainty)
whether there exists some ordering ``π`` of the rows and columns of ``A`` such that the
bandwidth of ``PAPᵀ`` is at most `k`, where ``P`` is the permutation matrix corresponding to
``π``. This is known to be decidable in ``O(nᵏ)`` time, although some deciders (e.g.,
[`CapraraSalazarGonzalez`](@ref)) run in exponential time instead to produce even quicker
runtimes in practice.

If ``k ≥ n - 1``, then this function immediately answers in the affirmative, since the
maximum possible bandwidth of an ``n×n`` matrix is ``n - 1``. After this initial check, a
preliminary lower bound on the bandwidth is computed in ``O(n³)`` time using results from
Caprara and Salazar-González (2005). If this lower bound is greater than ``k```

# Arguments
- `A::AbstractMatrix{<:Number}`: the (square) matrix whose bandwidth is tested.
- `k::Integer`: the threshold bandwidth against which to test.
- `decider::AbstractDecider`: the matrix bandwidth recognition algorithm to use; defaults to
    [`CapraraSalazarGonzalez`](@ref). (See the [`Recognition`](@ref) module documentation
    for a full list of supported deciders.)

# Returns
- `::RecognitionResult`: a struct containing the algorithm used, the original matrix `A`,
    the identified ordering of the rows and columns (if one exists), the threshold bandwidth
    `k`, and a boolean indicating whether the ordering exists.

# Examples
Multiple algorithms to decide whether a given matrix has bandwidth at most `k` are
available. Naturally, they will always agree, but the final orderings produced (in the case
of an affirmative) may differ:

```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(52);

julia> (n, p) = (8, 0.2);

julia> A = sprand(Bool, n, n, p);

julia> A = A .|| A' # Ensure structural symmetry
8×8 SparseMatrixCSC{Bool, Int64} with 22 stored entries:
 1  ⋅  ⋅  1  ⋅  1  1  ⋅
 ⋅  ⋅  ⋅  ⋅  1  1  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1
 1  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅
 ⋅  1  ⋅  ⋅  ⋅  ⋅  1  1
 1  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 1  ⋅  ⋅  1  1  ⋅  ⋅  1
 ⋅  ⋅  1  ⋅  1  ⋅  1  1

julia> k = 3;

julia> res_csg = has_bandwidth_k_ordering(A, k, Recognition.CapraraSalazarGonzalez())
Results of Bandwidth Recognition Algorithm
 * Algorithm: Caprara–Salazar-González
 * Bandwidth Threshold k: 3
 * Has Bandwidth ≤ k Ordering: true
 * Original Bandwidth: 6
 * Matrix Size: 8×8

julia> A[res_csg.ordering, res_csg.ordering]
8×8 SparseMatrixCSC{Bool, Int64} with 22 stored entries:
 ⋅  ⋅  1  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  1  ⋅  1  ⋅  ⋅  ⋅
 1  1  1  1  ⋅  ⋅  ⋅  ⋅
 1  ⋅  1  ⋅  ⋅  1  1  ⋅
 ⋅  1  ⋅  ⋅  ⋅  1  ⋅  ⋅
 ⋅  ⋅  ⋅  1  1  ⋅  1  ⋅
 ⋅  ⋅  ⋅  1  ⋅  1  1  1
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅

julia> res_sgs = has_bandwidth_k_ordering(A, k, Recognition.SaxeGurariSudborough())
Results of Bandwidth Recognition Algorithm
 * Algorithm: Saxe–Gurari–Sudborough
 * Bandwidth Threshold k: 3
 * Has Bandwidth ≤ k Ordering: true
 * Original Bandwidth: 6
 * Matrix Size: 8×8

julia> A[res_sgs.ordering, res_sgs.ordering]
8×8 SparseMatrixCSC{Bool, Int64} with 22 stored entries:
 ⋅  1  ⋅  1  ⋅  ⋅  ⋅  ⋅
 1  ⋅  1  ⋅  1  ⋅  ⋅  ⋅
 ⋅  1  1  ⋅  1  1  ⋅  ⋅
 1  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅
 ⋅  1  1  ⋅  ⋅  ⋅  1  1
 ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  1  1  ⋅  1  1
 ⋅  ⋅  ⋅  ⋅  1  ⋅  1  ⋅
```

If no decider is specified, then the Caprara–Salazar-González algorithm is used by default:

```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(174);

julia> (n, p, k) = (20, 0.1, 4);

julia> A = sprand(n, n, p);

julia> A = A .+ A' # Ensure structural symmetry;

julia> has_bandwidth_k_ordering(A, k)
Results of Bandwidth Recognition Algorithm
 * Algorithm: Caprara–Salazar-González
 * Bandwidth Threshold k: 4
 * Has Bandwidth ≤ k Ordering: false
 * Original Bandwidth: 15
 * Matrix Size: 20×20
```

# Notes
To implement a new matrix bandwidth recognition algorithm, define a new concrete subtype of
[`AbstractDecider`](@ref) then implement a corresponding
`_has_bandwidth_k_ordering_impl(::AbstractMatrix{Bool}, ::Integer, ::NewDeciderType)`
method. Do *not* attempt to directly implement a new `has_bandwidth_k_ordering` method, as
the function contains common preprocessing logic independent of the specific algorithm used.

Note also that some texts define matrix bandwidth to be the minimum non-negative integer
``k`` such that ``A[i, j] = 0`` whenever ``|i - j| ≥ k`` instead, particularly in more
mathematically-minded communities. Effectively, this definition treats diagonal matrices as
bandwidth ``1``, tridiagonal matrices as bandwidth ``2``, and so on. Our definition, on the
other hand, is more common in computer science contexts, treating diagonal matrices as
bandwidth ``0`` and tridiagonal matrices as bandwidth ``1``. (Both definitions, however,
agree that the bandwidth of an empty matrix is simply ``0``.)
"""
function has_bandwidth_k_ordering(
    A::AbstractMatrix{<:Number}, k::Integer, decider::AbstractDecider=DEFAULT_DECIDER
)
    if !allequal(size(A))
        throw(RectangularMatrixError(A))
    end

    if _requires_structural_symmetry(decider) && !is_structurally_symmetric(A)
        throw(StructuralAsymmetryError(A, decider))
    end

    #= If the bandwidth of `A` is already less than or equal to `k`, then the current
    ordering suffices. Otherwise, we compute a preliminary lower bound on the bandwidth
    using results from Caprara and Salazar-González (2005) to determine whether we should
    continue searching for an ordering. =#
    if bandwidth(A) <= k
        ordering = collect(axes(A, 1))
    elseif bandwidth_lower_bound(A) <= k
        #= We are only concerned with which (off-diagonal) entries are nonzero, not the actual
        values. We also set every diagonal entry to `false` for consistency with any algorithms
        that assume an adjacency matrix structure. =#
        A_bool = offdiag_nz_support(A)
        ordering = _has_bandwidth_k_ordering_impl(A_bool, k, decider)
    else
        ordering = nothing
    end

    return RecognitionResult(decider, A, ordering, k)
end

#= Compute an ordering inducing a bandwidth of at most `k` for a preprocessed
`AbstractMatrix{Bool}`, if one exists; otherwise, return `nothing`. Restricting entries to
booleans can improve performance via cache optimizations, bitwise operations, etc. Each
concrete subtype of `AbstractDecider` must implement its own
`_has_bandwidth_k_ordering_impl` method to define the corresponding algorithm logic. =#
function _has_bandwidth_k_ordering_impl(
    ::AbstractMatrix{Bool}, ::Integer, ::T
) where {T<:AbstractDecider}
    throw(NotImplementedError(_has_bandwidth_k_ordering_impl, :decider, T, AbstractDecider))
end
