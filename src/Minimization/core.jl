# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    minimize_bandwidth(A, solver=ReverseCuthillMcKee()) -> MinimizationResult

Minimize the bandwidth of `A` using the algorithm defined by `solver`.

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

This function computes a (near-)optimal ordering ``π`` of the rows and columns of ``A`` so
that the bandwidth of ``PAPᵀ`` is minimized, where ``P`` is the permutation matrix
corresponding to ``π``. This is known to be an NP-complete problem; however, several
heuristic algorithms such as [`ReverseCuthillMcKee`](@ref) run in polynomial time while
still producing near-optimal orderings in practice. Exact methods like
[`CapraraSalazarGonzalez`](@ref) are also available, but they are exponential in time
complexity and thus only feasible for relatively small matrices.

# Arguments
- `A::AbstractMatrix{<:Number}`: the (square) matrix whose bandwidth is minimized.
- `solver::AbstractSolver`: the matrix bandwidth minimization algorithm to use; defaults to
    [`ReverseCuthillMcKee`](@ref). (See the [`Minimization`](@ref) module documentation for
    a full list of supported solvers.)

# Returns
- `::MinimizationResult`: a struct containing the original matrix `A`, the minimized
    bandwidth, the (near-)optimal ordering of the rows and columns, and the algorithm used.

# Examples
[TODO: Add here once more solvers are implemented]

# Notes
Some texts define matrix bandwidth to be the minimum non-negative integer ``k`` such that
``A[i, j] = 0`` whenever ``|i - j| ≥ k`` instead, particularly in more mathematically-minded
communities. Effectively, this definition treats diagonal matrices as bandwidth ``1``,
tridiagonal matrices as bandwidth ``2``, and so on. Our definition, on the other hand, is
more common in computer science contexts, treating diagonal matrices as bandwidth ``0`` and
tridiagonal matrices as bandwidth ``1``. (Both definitions, however, agree that the
bandwidth of an empty matrix is simply ``0``.)
"""
function minimize_bandwidth(
    A::AbstractMatrix{<:Number}, solver::AbstractSolver=DEFAULT_SOLVER
)
    if !allequal(size(A))
        throw(RectangularMatrixError(A))
    end

    if _requires_symmetry(solver) && !_is_structurally_symmetric(A)
        throw(StructuralAsymmetryError(A, solver))
    end

    #= We are only concerned with which (off-diagonal) entries are nonzero, not the actual
    values. We also set every diagonal entry to `false` for consistency with any algorithms
    that assume an adjacency matrix structure. =#
    A_bool = _offdiag_nonzero_support(A)

    bandwidth_orig = bandwidth(A_bool)

    # If `A` is already diagonal/empty, no possible ordering can improve its bandwidth
    if bandwidth_orig == 0
        bandwidth_min = 0
        ordering = collect(axes(A_bool, 1)) # The original ordering
    else
        #= We call the `_bool_minimal_band_ordering` helper function, which is the function
        on which dispatch is actually defined for each concrete `AbstractSolver` subtype. =#
        ordering = _bool_minimal_band_ordering(A_bool, solver)
        A_reordered = A_bool[ordering, ordering]
        bandwidth_min = bandwidth(A_reordered)

        #= Some heuristic and metaheuristic algorithms may induce a bandwidth larger than
        the original if the input matrix is already optimally ordered (or close to it). =#
        if bandwidth_min >= bandwidth_orig
            bandwidth_min = bandwidth_orig
            ordering = collect(axes(A_bool, 1)) # The original ordering
        end
    end

    return MinimizationResult(solver, A, ordering, bandwidth_min)
end

#= Compute a minimal bandwidth ordering for a preprocessed `AbstractMatrix{Bool}`.
Restricting entries to booleans can improve performance via cache optimizations, bitwise
operations, etc. Each concrete subtype of `AbstractSolver` must implement its own
`_bool_minimal_band_ordering` method to define the corresponding algorithm logic. =#
function _bool_minimal_band_ordering(::AbstractMatrix{Bool}, ::T) where {T<:AbstractSolver}
    throw(NotImplementedError(_bool_minimal_band_ordering, :solver, T, AbstractSolver))
end
