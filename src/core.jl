# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    minimize_bandwidth(A, solver=ReverseCuthillMcKee()) -> BandwidthResult

Minimize the matrix bandwidth of `A` using the algorithm defined by `solver`.

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

This function computes a (near-)optimal ordering ``π`` of the rows and columns of ``A`` so
that the bandwidth of ``PAPᵀ`` is minimized, where ``P`` is the permutation matrix
corresponding to ``π``. This is known to be an NP-complete problem; however, several
heuristic algorithms such as [`ReverseCuthillMcKee`](@ref) run in polynomial time while
still producing near-optimal orderings in practice. Exact methods like [`MBPS`](@ref) are
also available, but they are exponential in time complexity and thus only feasible for
relatively small matrices.

# Arguments
- `A::AbstractMatrix{T}`: the (square) matrix whose bandwidth is to be minimized.
- `solver::AbstractSolver`: the matrix bandwidth minimization algorithm to use; defaults to
    [`ReverseCuthillMcKee`](@ref). (See the [`MatrixBandwidth`](@ref) module documentation
    for a full list of supported solvers.)

# Returns
- `::BandwidthResult`: a struct containing the original matrix `A`, the minimized bandwidth,
    the (near-)optimal ordering of the rows and columns, and the algorithm used.

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
    A::AbstractMatrix{T}, solver::AbstractSolver=DEFAULT_SOLVER
) where {T<:Number}
    _assert_matrix_is_square(A) # Bandwidth is not defined for non-square matrices

    #= Converting entries to booleans can improve performance via cache optimizations,
    bitwise operations, etc. (If `A` is already an `AbstractMatrix{Bool}`, we simply alias
    to `A_bool` rather than copying to avoid reallocation overhead.) =#
    if T === Bool
        A_bool = A
    else
        A_bool = (!iszero).(A)
    end

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

    return BandwidthResult(A, bandwidth_min, ordering, solver)
end

"""
    bandwidth(A) -> Int

Compute the bandwidth of `A` without any permutations.

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

In contrast to [`minimize_bandwidth`](@ref), this function does not attempt to minimize the
bandwidth of `A` by permuting its rows and columns—it simply computes its bandwidth as is.

# Arguments
- `A::AbstractMatrix{T}`: the (square) matrix whose bandwidth is to be computed.

# Returns
- `::Int`: the bandwidth of `A`.

# Examples
`bandwidth` correctly identifies the bandwidth of a pentadiagonal matrix as ``2`` and does
not attempt to find a minimizing permutation upon shuffling of its rows and columns:
```jldoctest
julia> using Random

julia> Random.seed!(242622);

julia> (n, k) = (8, 2);

julia> perm = randperm(n);

julia> A = (!iszero).(random_banded_matrix(8, 2))
8×8 BitMatrix:
 1  0  1  0  0  0  0  0
 1  0  1  1  0  0  0  0
 1  1  1  1  1  0  0  0
 0  1  0  1  1  1  0  0
 0  0  0  1  1  1  1  0
 0  0  0  1  1  1  1  1
 0  0  0  0  1  1  1  1
 0  0  0  0  0  1  0  1

julia> bandwidth(A)
2

julia> A_shuffled = A[perm, perm]
8×8 BitMatrix:
 1  1  1  0  0  0  1  1
 1  1  1  0  0  0  0  1
 1  1  1  0  0  0  1  0
 0  0  1  1  1  1  1  0
 0  0  0  1  1  0  0  0
 0  0  0  1  1  0  1  0
 1  0  1  0  0  1  1  0
 1  0  0  0  0  0  0  1

julia> bandwidth(A_shuffled)
7
```

# Notes
Some texts define matrix bandwidth to be the minimum non-negative integer ``k`` such that
``A[i, j] = 0`` whenever ``|i - j| ≥ k`` instead, particularly in more mathematically-minded
communities. Effectively, this definition treats diagonal matrices as bandwidth ``1``,
tridiagonal matrices as bandwidth ``2``, and so on. Our definition, on the other hand, is
more common in computer science contexts, treating diagonal matrices as bandwidth ``0`` and
tridiagonal matrices as bandwidth ``1``. (Both definitions, however, agree that the
bandwidth of an empty matrix is simply ``0``.)
"""
function bandwidth(A::AbstractMatrix{T}) where {T<:Number}
    _assert_matrix_is_square(A)

    if T === Bool
        A_bool = A
    else
        A_bool = (!iszero).(A)
    end

    indices = findall(A_bool)

    if isempty(indices)
        band = 0
    else
        band = maximum(abs(index[1] - index[2]) for index in indices)
    end

    return band
end

#= Compute a minimal bandwidth ordering for a preprocessed `AbstractMatrix{Bool}`.
Restricting entries to booleans can improve performance via cache optimizations, bitwise
operations, etc. Each concrete `AbstractSolver` subtype must implement its own
`_bool_minimal_band_ordering` method to define the corresponding algorithm logic. =#
function _bool_minimal_band_ordering(::AbstractMatrix{Bool}, ::T) where {T<:AbstractSolver}
    S = supertype(T)

    if S === AbstractSolver
        subtype = T
    else
        subtype = S
    end

    throw(
        NotImplementedError(_bool_minimal_band_ordering, :solver, subtype, AbstractSolver)
    )
end
