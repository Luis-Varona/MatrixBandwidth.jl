# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    minimize_bandwidth(A, solver=GibbsPooleStockmeyer()) -> MinimizationResult

Minimize the bandwidth of `A` using the algorithm defined by `solver`.

The *bandwidth* of an ``n×n`` matrix ``A`` is the minimum non-negative integer
``k ∈ \\{0, 1, …, n - 1\\}`` such that ``A[i, j] = 0`` whenever ``|i - j| > k``.

This function computes a (near-)optimal ordering ``π`` of the rows and columns of ``A`` so
that the bandwidth of ``PAPᵀ`` is minimized, where ``P`` is the permutation matrix
corresponding to ``π``. This is known to be an NP-complete problem; however, several
heuristic algorithms such as Gibbs–Poole–Stockmeyer run in polynomial time while still
still producing near-optimal orderings in practice. Exact methods like
Caprara–Salazar-González are also available, but they are at least exponential in time
complexity and thus only feasible for relatively small matrices.

# Arguments
- `A::AbstractMatrix{<:Number}`: the (square) matrix whose bandwidth is minimized.
- `solver::AbstractSolver`: the matrix bandwidth minimization algorithm to use; defaults to
    [`GibbsPooleStockmeyer`](@ref). (See the [`Minimization`](@ref) module documentation for
    a full list of supported solvers.)

# Returns
- `::MinimizationResult`: a struct containing the algorithm used, the original matrix `A`,
    the (near-)optimal ordering of the rows and columns, and the minimized bandwidth.

# Examples
Multiple algorithms to minimize the bandwidth of a given matrix are available. In
particular, there are exact solvers (which always guarantee optimal solutions), heuristic
solvers (which produce near-optimal solutions more quickly than exact methods), and
metaheuristic solvers (which employ iterative search frameworks to find better solutions
than heuristic methods, but even more slowly).

Certainly, exact solvers will always produce the same optimal bandwidth (but likely
different orderings):

```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(38);

julia> (n, p) = (20, 0.05);

julia> A = sprand(n, n, p);

julia> A = A .+ A' # Ensure structural symmetry
20×20 SparseMatrixCSC{Float64, Int64} with 36 stored entries:
⎡⠀⠀⠀⠀⠀⣎⠁⠈⠀⠂⎤
⎢⠀⠀⠀⠀⡀⠐⠀⠀⠀⠠⎥
⎢⡠⢤⢀⠈⠠⠂⠀⠀⠔⠀⎥
⎢⡁⠀⠀⠀⠀⠀⠀⠀⡀⠐⎥
⎣⠠⠀⠀⡀⠐⠁⢀⠈⢐⠔⎦

julia> res_dcm = minimize_bandwidth(A, Minimization.DelCorsoManzini())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Del Corso–Manzini
 * Approach: exact
 * Minimum Bandwidth: 3
 * Original Bandwidth: 17
 * Matrix Size: 20×20

julia> A[res_dcm.ordering, res_dcm.ordering]
20×20 SparseMatrixCSC{Float64, Int64} with 36 stored entries:
⎡⠠⠂⠤⡀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠣⡁⠈⠀⡀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠠⡤⠋⠤⡀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠣⡀⡨⠢⡀⎥
⎣⠀⠀⠀⠀⠀⠀⠈⠢⠁⠀⎦

julia> res_sgs = minimize_bandwidth(A, Minimization.SaxeGurariSudborough())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Saxe–Gurari–Sudborough
 * Approach: exact
 * Minimum Bandwidth: 3
 * Original Bandwidth: 17
 * Matrix Size: 20×20

julia> A[res_sgs.ordering, res_sgs.ordering]
20×20 SparseMatrixCSC{Float64, Int64} with 36 stored entries:
⎡⠀⡠⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠔⣡⠢⡀⠀⠀⠀⠀⎥
⎢⠀⠀⠈⠢⡀⡨⣀⡀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠸⡀⠈⢂⡀⎥
⎣⠀⠀⠀⠀⠀⠀⠈⠰⠄⡡⎦
```

However, the answers of (meta)heuristic solvers may differ from each other:

```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(405);

julia> (n, p) = (400, 0.002);

julia> A = sprand(n, n, p);

julia> A = A .+ A' # Ensure structural symmetry;

julia> minimize_bandwidth(A, Minimization.GibbsPooleStockmeyer())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Gibbs–Poole–Stockmeyer
 * Approach: heuristic
 * Minimum Bandwidth: 31
 * Original Bandwidth: 380
 * Matrix Size: 400×400

julia> minimize_bandwidth(A, Minimization.ReverseCuthillMcKee())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Reverse Cuthill–McKee
 * Approach: heuristic
 * Minimum Bandwidth: 37
 * Original Bandwidth: 380
 * Matrix Size: 400×400
```

If no solver is specified, then the heuristic Gibbs–Poole–Stockmeyer algorithm is used by
default:

```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(80);

julia> (n, p) = (500, 0.001);

julia> A = sprand(n, n, p);

julia> A = A .+ A' # Ensure structural symmetry
500×500 SparseMatrixCSC{Float64, Int64} with 496 stored entries:
⎡⠀⠀⠁⠐⠀⠀⠀⠄⠠⢐⠄⠀⡈⠀⠀⡆⠠⠑⠀⠰⠀⠀⠀⠐⠀⠄⠀⠀⡠⠀⠀⠀⠆⠐⠀⠄⢠⡈⠈⠀⎤
⎢⢁⠀⠄⠁⠠⠊⠀⠀⠂⠀⡀⠂⠀⢀⠀⠀⠀⠈⠠⠠⠀⠠⠠⠀⠀⡀⠀⠁⡀⠀⡀⠐⠀⠄⠀⠀⠄⠂⡈⠠⎥
⎢⠀⠀⡠⠂⠅⠁⠀⠀⠀⠀⠐⡁⠀⡀⠈⢈⠀⠐⠀⠀⠀⠀⠃⢀⠀⠀⠀⠀⠀⠀⠐⠀⠤⠁⠀⠁⠀⠀⠀⡀⎥
⎢⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠀⠠⠀⠂⠁⠀⠀⡐⠔⠂⠂⠀⠠⠀⠐⠀⠀⠄⠀⠄⠀⠀⡀⠀⢀⠁⠀⎥
⎢⢀⢂⠈⠀⠀⠀⠀⠀⠎⠁⠄⠂⠀⠁⣀⠁⠈⠠⠀⠀⠂⠄⢀⠄⠀⠠⠀⠀⠉⠈⠀⠀⠀⠀⠠⠀⠀⠀⠀⡀⎥
⎢⠀⠁⠠⠈⠔⠠⠠⠀⠠⠁⠊⠄⠀⢀⠀⠀⠠⠀⠠⡌⠀⠈⠀⠀⠀⠠⠈⠀⠀⠀⠂⠀⠀⠀⠀⠡⠀⠠⠠⢂⎥
⎢⠂⠈⠀⢀⠀⠠⠀⠀⠄⠀⠀⢀⡐⠈⠀⠀⠀⠀⠀⠈⠐⠠⠀⠀⠀⠀⠀⡠⠑⠀⠂⠀⠀⠀⡀⠀⡐⠀⠄⠁⎥
⎢⠠⠤⠀⠀⡂⢀⠀⠂⠄⠘⠀⠀⠀⠀⠀⠀⠀⠑⠀⠂⡄⠀⠄⠀⡁⡐⠁⠠⠂⠀⠀⠐⠀⠀⠁⠀⠀⠁⠈⠀⎥
⎢⢄⠂⡀⠀⢀⠀⠌⠀⠂⡀⠀⠂⠀⠀⢄⠀⠀⠀⠀⢀⠀⠀⢀⠄⡀⠂⠂⠀⠅⠀⠀⠀⠀⠁⠀⢀⠄⠀⠀⠀⎥
⎢⢀⡀⠀⡂⠀⠀⠀⠀⠀⠀⡀⠦⡀⠀⠠⠀⠀⢀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⡐⢀⢀⠁⠔⠆⠀⠀⠀⢀⠅⎥
⎢⠀⠀⠀⡀⠀⠀⢐⠌⠈⠄⡀⠀⠐⡀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀⠀⠌⡀⠀⠀⠀⠀⠀⠀⠂⎥
⎢⢀⠀⠀⠂⠉⢀⠨⠀⠀⠔⠀⠀⠀⠀⠀⠁⠀⠔⠐⠀⠀⠀⠀⡠⠠⠀⠰⠀⠠⠘⠀⠐⠀⡀⠁⠠⠀⡑⠀⠀⎥
⎢⠀⠄⠀⠠⠀⠀⠀⡀⠀⡀⠀⡀⠀⠀⢁⠨⠠⠈⠀⠀⠀⠀⠀⠂⠀⠀⠁⠀⠁⠂⠄⠀⠠⠀⠀⠀⠀⠀⡀⠀⎥
⎢⠀⠀⠄⠀⠀⠀⢀⠀⠀⠀⠂⠀⠀⡠⠁⡀⠈⠀⠀⠀⠀⠠⠐⠂⠁⠀⠀⠀⢁⠊⠠⠀⢀⠀⠀⠀⠰⠌⠀⠀⎥
⎢⠀⠊⠀⠈⠀⠀⠀⠀⡃⠀⠀⠀⠑⠀⠈⠀⠁⠁⢀⠠⠀⠀⣀⠂⠡⠀⡡⠐⠀⠀⠀⠐⠂⠰⠀⠠⠀⠀⠂⠀⎥
⎢⠀⠀⢀⠈⠐⠀⠀⠁⠀⠀⠈⠀⠈⠀⢀⠀⠀⠀⠀⢐⡀⠄⢀⠀⠀⠁⠀⠂⢀⠀⠤⠃⠀⠐⠀⠀⠀⠈⠀⠀⎥
⎢⢈⠁⠀⠄⠄⠃⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⢁⠄⠀⠈⠀⠠⠀⠂⠀⠐⢈⡀⢀⠀⠄⠁⠀⠈⠀⠆⠀⠀⎥
⎢⠀⠄⠀⠀⠄⠀⠀⠠⠀⠂⠄⡀⠀⠈⠁⠀⠀⢀⠈⠁⠀⠀⠁⡀⠀⠀⠀⠀⠀⡀⠀⠀⡀⠀⢠⠒⠈⠀⠀⠀⎥
⎢⡀⠲⠠⠁⠀⠀⠀⢀⠀⠀⠀⡀⠐⠈⠄⠀⠀⠁⠀⠀⠀⠀⢄⠠⠀⠀⡐⠆⠀⠀⡀⠀⠠⠄⠂⠀⠤⠃⠀⠁⎥
⎣⠂⠀⠂⡈⠀⠠⠁⠀⠀⠠⠠⢂⠄⠁⠂⠀⠀⠀⠄⠔⠠⠀⠀⠀⠀⠈⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠄⠀⠁⠀⎦

julia> res = minimize_bandwidth(A)
Results of Bandwidth Minimization Algorithm
 * Algorithm: Gibbs–Poole–Stockmeyer
 * Approach: heuristic
 * Minimum Bandwidth: 6
 * Original Bandwidth: 487
 * Matrix Size: 500×500

julia> A[res.ordering, res.ordering]
500×500 SparseMatrixCSC{Float64, Int64} with 496 stored entries:
⎡⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠑⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠐⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠊⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠻⣢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠱⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠱⢆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠱⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⢆⡀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⢆⡀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⎦
```

# Notes
To implement a new matrix bandwidth minimization algorithm, define a new concrete subtype of
[`AbstractSolver`](@ref) (or of one of its abstract subtypes like
[`MetaheuristicSolver`](@ref)) then implement a corresponding
`_minimize_bandwidth_impl(::AbstractMatrix{Bool}, ::NewSolverType)` method. Do *not* attempt
to directly implement a new `minimize_bandwidth` method, as the function contains common
preprocessing logic independent of the specific algorithm used.

Note also that some texts define matrix bandwidth to be the minimum non-negative integer
``k`` such that ``A[i, j] = 0`` whenever ``|i - j| ≥ k`` instead, particularly in more
mathematically-minded communities. Effectively, this definition treats diagonal matrices as
bandwidth ``1``, tridiagonal matrices as bandwidth ``2``, and so on. Our definition, on the
other hand, is more common in computer science contexts, treating diagonal matrices as
bandwidth ``0`` and tridiagonal matrices as bandwidth ``1``. (Both definitions, however,
agree that the bandwidth of an empty matrix is simply ``0``.)
"""
function minimize_bandwidth(
    A::AbstractMatrix{<:Number}, solver::AbstractSolver=DEFAULT_SOLVER
)
    if !allequal(size(A))
        throw(RectangularMatrixError(A))
    end

    if _requires_structural_symmetry(solver) && !is_structurally_symmetric(A)
        throw(StructuralAsymmetryError(A, solver))
    end

    #= We are only concerned with which (off-diagonal) entries are nonzero, not the actual
    values. We also set every diagonal entry to `false` for consistency with any algorithms
    that assume an adjacency matrix structure. =#
    A_bool = offdiag_nz_support(A)

    bandwidth_orig = bandwidth(A_bool)

    # If `A` is already optimally ordered, no ordering can further reduce its bandwidth
    if bandwidth_orig == bandwidth_lower_bound(A_bool)
        bandwidth_min = bandwidth_orig
        ordering = collect(axes(A_bool, 1)) # The original ordering
    else
        #= We call the `_minimize_bandwidth_impl` helper function, which is the function on
        which dispatch is actually defined for each concrete `AbstractSolver` subtype. =#
        ordering = _minimize_bandwidth_impl(A_bool, solver)
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
`_minimize_bandwidth_impl` method to define the corresponding algorithm logic. =#
function _minimize_bandwidth_impl(::AbstractMatrix{Bool}, ::T) where {T<:AbstractSolver}
    throw(NotImplementedError(_minimize_bandwidth_impl, :solver, T, AbstractSolver))
end
