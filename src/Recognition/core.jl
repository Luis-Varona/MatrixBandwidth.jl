# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    has_bandwidth_k_ordering(A, k, decider=CapraraSalazarGonzalez()) -> RecognitionResult

Determine whether `A` has bandwidth at most `k` using the algorithm defined by `decider`.

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

This function [TODO: Write here]

# Arguments
- `A::AbstractMatrix{<:Number}`: the (square) matrix whose bandwidth is tested.
- `k::Int`: the threshold bandwidth against which to test.
- `decider::AbstractDecider`: the matrix bandwidth recognition algorithm to use; defaults to
    [`CapraraSalazarGonzalez`](@ref). (See the [`Recognition`](@ref) module documentation
    for a full list of supported deciders.)

# Returns
- `::RecognitionResult`: TODO: Write here

# Examples
[TODO: Add here once more deciders are implemented]

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
    A::AbstractMatrix{<:Number}, k::Int, decider::AbstractDecider=DEFAULT_DECIDER
)
    if !allequal(size(A))
        throw(RectangularMatrixError(A))
    end

    if _requires_symmetry(decider) && !_is_structurally_symmetric(A)
        throw(StructuralAsymmetryError(A, decider))
    end

    #= We are only concerned with which (off-diagonal) entries are nonzero, not the actual
    values. We also set every diagonal entry to `false` for consistency with any algorithms
    that assume an adjacency matrix structure. =#
    A_bool = _offdiag_nonzero_support(A)

    bandwidth_orig = bandwidth(A_bool)

    #= If the bandwidth of `A` is already less than or equal to `k`, then the current
    row/column ordering suffices. =#
    if bandwidth_orig <= k
        bandwidth_k_ordering = collect(axes(A_bool, 1)) # The original ordering
    else
        #= Compute a preliminary lower bound on the bandwidth using results from Caprara and
        Salazar-González (2005). =#
        lower_bound = bandwidth_lower_bound(A_bool)

        if lower_bound > k
            bandwidth_k_ordering = nothing
        else # The default case wherein a more expensive exact decider is used
            bandwidth_k_ordering = _bool_bandwidth_k_ordering(A_bool, k, decider)
        end
    end

    return RecognitionResult(decider, A, bandwidth_k_ordering, k)
end

# TODO: Summarize here. Returns either `nothing` or a `Vector{Int}`.
function _bool_bandwidth_k_ordering(
    ::AbstractMatrix{Bool}, ::Int, ::T
) where {T<:AbstractDecider}
    throw(NotImplementedError(_bool_bandwidth_k_ordering, :decider, T, AbstractDecider))
end
