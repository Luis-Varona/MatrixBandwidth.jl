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
Equivalently, ``A`` has bandwidth *at most* ``k`` if all entries above the ``k``ᵗʰ
superdiagonal and below the ``k``ᵗʰ subdiagonal are zero, and ``A`` has bandwidth *at least*
``k`` if there exists any nonzero entry in the ``k``ᵗʰ superdiagonal or subdiagonal.

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
[TODO: Add here once more deciders are implemented. For now, refer to the **Examples**
sections of the [`DelCorsoManzini`](@ref) and [`DelCorsoManziniWithPS`](@ref) docstrings.]

# Notes
Some texts define matrix bandwidth to be the minimum non-negative integer ``k`` such that
``A[i, j] = 0`` whenever ``|i - j| ≥ k`` instead, particularly in more mathematically-minded
communities. Effectively, this definition treats diagonal matrices as bandwidth ``1``,
tridiagonal matrices as bandwidth ``2``, and so on. Our definition, on the other hand, is
more common in computer science contexts, treating diagonal matrices as bandwidth ``0`` and
tridiagonal matrices as bandwidth ``1``. (Both definitions, however, agree that the
bandwidth of an empty matrix is simply ``0``.)
"""
function has_bandwidth_k_ordering(
    A::AbstractMatrix{<:Number}, k::Integer, decider::AbstractDecider=DEFAULT_DECIDER
)
    if !allequal(size(A))
        throw(RectangularMatrixError(A))
    end

    if _requires_structural_symmetry(decider) && !_is_structurally_symmetric(A)
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
        A_bool = _offdiag_nonzero_support(A)
        ordering = _bool_bandwidth_k_ordering(A_bool, k, decider)
    else
        ordering = nothing
    end

    return RecognitionResult(decider, A, ordering, k)
end

#= Compute an ordering inducing a bandwidth of at most `k` for a preprocessed
`AbstractMatrix{Bool}`, if one exists; otherwise, return `nothing`. Restricting entries to
booleans can improve performance via cache optimizations, bitwise operations, etc. Each
concrete subtype of `AbstractDecider` must implement its own `_bool_bandwidth_k_ordering`
method to define the corresponding algorithm logic. =#
function _bool_bandwidth_k_ordering(
    ::AbstractMatrix{Bool}, ::Integer, ::T
) where {T<:AbstractDecider}
    throw(NotImplementedError(_bool_bandwidth_k_ordering, :decider, T, AbstractDecider))
end
