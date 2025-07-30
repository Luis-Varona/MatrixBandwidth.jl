# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CuthillMcKee <: HeuristicSolver <: AbstractSolver <: AbstractAlgorithm

The *Cuthill–McKee algorithm* is a heuristic method for minimizing the bandwidth of a
structurally symmetric matrix ``A``. It considers the graph ``G(A)`` whose adjacency matrix
is ``A`` (ignoring weights and self-loops) and performs a breadth-first search of each
connected component of ``G(A)``, starting from a low-degree node then visiting its neighbors
in order of increasing degree. Particularly effective when ``A`` is sparse, this heuristic
typically produces an ordering which induces a matrix bandwidth either equal to or very
close to the true minimum [CM69; pp. 157--58](@cite).

As noted above, the Cuthill–McKee algorithm requires structurally symmetric input (that is,
``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for ``1 ≤ i, j ≤ n``).

# Fields
- `node_selector::Function`: a function that selects a node from some connected component of
    the input matrix from which to start the breadth-first search. If no custom heuristic is
    specified, this field defaults to [`pseudo_peripheral_node`](@ref), which picks a node
    "farthest" from the others in the component (not necessarily the lowest-degree node).

# Supertype Hierarchy
`CuthillMcKee` <: [`HeuristicSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix ``A``, the Cuthill–McKee algorithm runs in ``O(n²)`` time.

[CG80](@cite) provide a linear-time implementation in the number of nonzero entries of
``A``, which is still quadratic when ``A`` is dense but often much faster when dealing with
sparse matrices. However, this would require that ``A`` be stored as a graph or a sparse
matrix, which runs counter to our desire to provide a bandwidth minimization API for all
`AbstractMatrix{<:Number}` types, including dense matrices. (In the future, however, we may
indeed consider supporting this more performant implementation for sparse matrices.)

It was found in [Geo71; pp. 114--15](@cite) that reversing the ordering produced by
Cuthill–McKee tends to induce a more optimal *matrix profile* (a measure of how far, on
average, nonzero entries are from the diagonal; see also [`MatrixBandwidth.profile`](@ref)).
This so-called *reverse Cuthill–McKee* variant is preferred in almost all cases—see
[`ReverseCuthillMcKee`](@ref) and the associated method of `_bool_minimal_band_ordering` for
    our implementation.

# Examples
In the following examples, [`MatrixBandwidth.random_banded_matrix`](@ref) is used to
generate random matrices with minimum bandwidth *close to* ``k``. In some cases, however,
the true minimum bandwidth up to symmetric permutation may be even less than ``k``, making
it hard to verify whether Cuthill–McKee finds a truly optimal ordering or simply a
near-optimal one. Nevertheless, the results are still very good in practice.

Cuthill–McKee finds a good ordering for a ``30×30`` matrix:
```jldoctest
julia> using Random

julia> Random.seed!(13);

julia> (n, k) = (30, 5);

julia> A = random_banded_matrix(n, k);

julia> perm = randperm(n);

julia> A_shuffled = A[perm, perm];

julia> bandwidth(A)
5

julia> bandwidth(A_shuffled) # Much larger after shuffling
25

julia> minimize_bandwidth(A_shuffled, Minimization.CuthillMcKee())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Cuthill–McKee
 * Approach: heuristic
 * Minimum Bandwidth: 5
 * Original Bandwidth: 25
 * Matrix Size: 30×30
```

Cuthill–McKee finds a good ordering for a structurally symmetric ``183×183`` matrix with
multiple (separate) connected components:
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(37452);

julia> (max_cc_size, max_band, p, num_ccs) = (60, 9, 0.2, 7);

julia> components = Vector{SparseMatrixCSC{Float64, Int64}}(undef, num_ccs);

julia> for i in 1:num_ccs # Some components may themselves be disconnected
           cc_size = rand(1:max_cc_size);
           cc_band = rand(0:min(max_band, cc_size - 1));
           components[i] = sparse(random_banded_matrix(cc_size, cc_band; p=p));
       end

julia> A = blockdiag(components...); # `A` has least 7 connected components

julia> perm = randperm(sum(map(cc -> size(cc, 1), components)));

julia> A_shuffled = A[perm, perm];

julia> res = minimize_bandwidth(A_shuffled, Minimization.CuthillMcKee());

julia> A # The original matrix
276×276 SparseMatrixCSC{Float64, Int64} with 464 stored entries:
⎡⢾⡷⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠘⢻⣲⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠘⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠈⠿⡧⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠉⢯⡷⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠚⣤⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢻⣶⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠯⡧⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⣤⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠛⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠱⣢⡀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⡢⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢴⣷⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠿⣧⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢿⡷⎦

julia> A_shuffled # A far-from-optimal ordering of `A`
276×276 SparseMatrixCSC{Float64, Int64} with 464 stored entries:
⎡⠁⢄⠀⢀⠀⠀⠀⢀⠠⠀⠀⠐⠀⠀⠀⠐⢀⡐⠀⠀⠀⢀⠀⠀⠀⠀⠐⠀⢠⠀⠀⠀⡄⠀⠀⠐⠀⠀⠂⠄⎤
⎢⠀⢀⠱⠂⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⢨⠀⠀⠀⠀⠀⡀⠁⠠⠀⠘⠀⠀⠡⢀⠈⠀⠀⠀⠀⠀⠀⠄⠀⠁⠁⎥
⎢⠀⠀⠀⠀⠑⢀⠀⠂⠀⠀⠀⠀⢐⠀⠀⠠⠈⠠⠀⠀⠀⠐⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⢢⢀⢀⠀⎥
⎢⠀⢀⡀⠀⠠⠀⠁⠄⠀⠠⠀⠄⠀⠀⠀⠄⠀⠀⠀⠀⢀⠀⠀⢀⠀⠑⠀⠀⠐⠠⠀⠀⠠⠨⠂⠀⠀⠀⠀⠀⎥
⎢⠀⠂⠀⠀⠀⠀⠀⡀⠱⢆⡀⠂⠀⠀⠀⠀⠀⠀⢀⢊⠀⠐⠐⠈⠀⠈⠀⢀⠄⠀⡀⠀⢁⢀⠠⠀⠃⠀⠊⠀⎥
⎢⢀⠀⠀⠀⠀⠀⠀⠄⠠⠈⠑⠀⢀⠐⠀⠌⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀⠀⢀⠉⢀⠀⠠⠈⠀⠀⣁⠁⎥
⎢⠀⠀⠀⠀⠐⠐⠀⠀⠀⠀⢀⠐⠁⠄⠈⠀⢌⠀⠆⠠⢀⠀⠄⠐⠰⠀⠀⠀⠁⠰⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⢀⠀⠂⠒⠀⡀⠀⠄⠀⠀⡀⠄⠂⠀⠐⢄⠁⢀⠀⠀⠀⡀⠀⠀⠀⠀⡠⠀⠀⠀⠀⠀⠀⠀⠀⢈⠀⠀⠀⠁⎥
⎢⢀⠰⠀⠀⠂⡀⠀⠀⠀⠀⠀⠈⠂⠑⠁⢀⠐⠄⠄⠂⠂⠜⠄⠀⠀⠀⡄⠀⠀⢀⠀⠠⠀⢀⠄⠀⢀⠀⠂⡂⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⡠⢐⠀⠀⠈⡁⠀⠀⠠⠁⠀⢀⠀⠀⠀⠀⡀⠀⠀⢀⠀⠈⠃⠀⠸⠈⠠⠀⠀⠀⢄⠂⎥
⎢⠀⢀⠄⠈⢀⠀⠀⠐⢀⠀⠀⠀⠀⠐⠀⠠⣈⠄⠀⠀⠐⢀⠀⡀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠊⠀⠠⠀⠐⠀⠀⎥
⎢⠀⠀⠀⠂⢀⠀⠀⢀⡐⠀⠀⠀⢀⠁⠀⠀⠀⠁⠀⠀⠀⠠⠄⣥⠉⠈⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⡀⠀⠀⠀⎥
⎢⠀⠀⠒⠀⠀⠀⢄⠀⡀⠀⠀⠀⠐⠂⠀⠀⠀⠀⠀⠈⠀⠀⡃⠀⠁⢀⠀⠀⢀⡀⢈⠈⠀⠀⠀⠂⠀⠠⠂⠂⎥
⎢⠐⠀⠄⡀⠀⠀⠀⠀⠀⢀⠀⠈⠀⠀⠀⠊⠀⠉⠀⢀⠀⠀⠀⠀⠀⠀⠑⠀⠀⠀⠀⠀⢀⠀⠈⠀⠛⠃⢄⠀⎥
⎢⠀⠒⡀⠐⠀⠀⠐⡀⠀⠁⠀⠀⢁⡀⠀⠀⠀⢀⡀⠀⠀⠀⠀⠀⠀⠰⠀⠀⠀⢄⠀⠰⠀⠠⠠⢀⠀⠀⢂⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⡄⠐⠂⠀⠀⠀⠀⡀⠉⠀⠐⠀⠀⠈⡂⠐⠀⠀⢀⡀⠀⣠⠀⠄⠠⠀⠀⡀⠀⠀⎥
⎢⠀⠉⠀⠀⠀⠀⡀⡂⠁⢐⠀⠐⠀⠀⠀⠀⠀⢀⡒⠂⡠⠀⠀⠀⠀⠀⠀⠐⠀⡀⠀⠄⠑⠄⠀⠀⠀⠀⠀⠀⎥
⎢⢀⠀⠀⠀⠀⡀⠈⠀⠀⠂⡀⠂⠀⠀⡀⢀⠀⠁⠀⠂⠀⡀⠀⠀⠠⠀⠂⠀⠀⢂⠀⠂⠀⠀⠁⢀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠁⠈⢒⠀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⢀⠀⠀⠈⠀⡀⠿⠀⠀⠀⠀⠠⠀⠀⠀⠀⠱⠆⠀⠀⎥
⎣⠈⠄⠅⠀⠀⠐⠀⠀⠊⠀⠅⠘⠀⠀⠄⠀⠨⠠⠠⠑⠀⠀⠀⠀⠨⠀⠀⠑⠈⠐⠀⠀⠀⠀⠀⠀⠀⠀⠔⢅⎦

julia> A_shuffled[res.ordering, res.ordering] # A near-optimal reordering of `A_shuffled`
276×276 SparseMatrixCSC{Float64, Int64} with 464 stored entries:
⎡⠱⣦⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠉⠻⣦⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠘⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⡦⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⡦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠺⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠚⣤⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⣤⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢄⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢄⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⎦

julia> bandwidth(A)
7

julia> bandwidth(A_shuffled) # Much larger after shuffling
266

julia> res # Even better than the original bandwidth (which was, clearly, not yet optimal)
Results of Bandwidth Minimization Algorithm
 * Algorithm: Cuthill–McKee
 * Approach: heuristic
 * Minimum Bandwidth: 5
 * Original Bandwidth: 266
 * Matrix Size: 276×276
```

# Notes
Note that the `node_selector` field must be of the form
`(A::AbstractMatrix{Bool}) -> Integer` (i.e., it must take in an boolean matrix and return
an integer). If this is not the case, an `ArgumentError` is thrown upon construction.
"""
struct CuthillMcKee <: HeuristicSolver
    node_selector::Function

    function CuthillMcKee(node_selector::Function=DEFAULT_SELECTOR)
        _assert_valid_node_selector(node_selector)
        return new(node_selector)
    end
end

Base.summary(::CuthillMcKee) = "Cuthill–McKee"

_requires_symmetry(::CuthillMcKee) = true

"""
    ReverseCuthillMcKee <: HeuristicSolver <: AbstractSolver <: AbstractAlgorithm

The *reverse Cuthill–McKee algorithm* is a variant of the *Cuthill–McKee algorithm*—a
heuristic method for minimizing the bandwidth of a structurally symmetric matrix ``A``.
Cuthill–McKee considers the graph ``G(A)`` whose adjacency matrix is ``A`` (ignoring weights
and self-loops) and performs a breadth-first search of each connected component of ``G(A)``,
starting from a low-degree node then visiting its neighbors in order of increasing degree.
Particularly effective when ``A`` is sparse, this heuristic typically produces an ordering
which induces a matrix bandwidth either equal to or very close to the true minimum
[CM69; pp. 157--58](@cite). The reverse Cuthill–McKee algorithm simply reverses the ordering
produced by application of Cuthill–McKee; it was found in [Geo71; pp. 114--15](@cite) that
although the bandwidth remains the same, this tends to produce a more optimal *matrix
profile* (a measure of how far, on average, nonzero entries are from the diagonal; see also
[`MatrixBandwidth.profile`](@ref)).

As noted above, the reverse Cuthill–McKee algorithm requires structurally symmetric input
(that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for
``1 ≤ i, j ≤ n``).

# Performance
Given an ``n×n`` input matrix ``A``, the reverse Cuthill–McKee algorithm runs in ``O(n²)``
time.

[CG80](@cite) provide a linear-time implementation in the number of nonzero entries of
``A``, which is still quadratic when ``A`` is dense but often much faster when dealing with
sparse matrices. However, this would require that ``A`` be stored as a graph or a sparse
matrix, which runs counter to our desire to provide a bandwidth minimization API for all
`AbstractMatrix{<:Number}` types, including dense matrices. (In the future, however, we may
indeed consider supporting this more performant implementation for sparse matrices.)

# Fields
- `node_selector::Function`: a function that selects a node from some connected component of
    the input matrix from which to start the breadth-first search. If no custom heuristic is
    specified, this field defaults to [`pseudo_peripheral_node`](@ref), which picks a node
    "farthest" from the others in the component (not necessarily the lowest-degree node).

# Supertype Hierarchy
`ReverseCuthillMcKee` <: [`HeuristicSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

# Examples
In the following examples, [`MatrixBandwidth.random_banded_matrix`](@ref) is used to
generate random matrices with minimum bandwidth *close to* ``k``. In some cases, however,
the true minimum bandwidth up to symmetric permutation may be even less than ``k``, making
it hard to verify whether reverse Cuthill–McKee finds a truly optimal ordering or simply a
near-optimal one. Nevertheless, the results are still very good in practice.

Reverse Cuthill–McKee finds a good ordering for a ``35×35`` matrix:
```jldoctest
julia> using Random

julia> Random.seed!(87);

julia> (n, k) = (35, 3);

julia> A = random_banded_matrix(n, k);

julia> perm = randperm(n);

julia> A_shuffled = A[perm, perm];

julia> bandwidth(A)
3

julia> bandwidth(A_shuffled) # Much larger after shuffling
30

julia> minimize_bandwidth(A_shuffled, Minimization.ReverseCuthillMcKee())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Reverse Cuthill–McKee
 * Approach: heuristic
 * Minimum Bandwidth: 3
 * Original Bandwidth: 30
 * Matrix Size: 35×35
```

Reverse Cuthill–McKee finds a good ordering for a ``235×235`` matrix with multiple
(separate) connected components:
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(5747);

julia> (max_cc_size, max_band, p, num_ccs) = (60, 9, 0.2, 8);

julia> components = Vector{SparseMatrixCSC{Float64, Int64}}(undef, num_ccs);

julia> for i in 1:num_ccs # Some components may themselves be disconnected
           cc_size = rand(0:max_cc_size);
           cc_band = rand(1:min(max_band, cc_size - 1));
           components[i] = sparse(random_banded_matrix(cc_size, cc_band; p=p));
       end

julia> A = blockdiag(components...); # `A` has least 8 connected components

julia> perm = randperm(sum(map(cc -> size(cc, 1), components)));

julia> A_shuffled = A[perm, perm];

julia> res = minimize_bandwidth(A_shuffled, Minimization.ReverseCuthillMcKee());

julia> A # The original matrix
235×235 SparseMatrixCSC{Float64, Int64} with 445 stored entries:
⎡⢾⣳⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠘⢿⡷⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠏⣥⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠉⢴⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠙⠻⢂⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⣮⣿⣢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠚⢿⡳⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢰⣶⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠾⡧⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⣢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⡠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢻⠖⣀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢏⡱⣄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢮⣷⣄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢻⣲⣄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢻⢖⎦

julia> A_shuffled # A far-from-optimal ordering of `A`
235×235 SparseMatrixCSC{Float64, Int64} with 445 stored entries:
⎡⠑⠄⠀⠀⠀⠀⠀⢀⢀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⡐⠀⠐⠀⠂⠀⠀⠀⠀⠀⠀⡂⠀⢀⠄⠁⠠⠐⠀⎤
⎢⠀⠀⠀⢄⡀⠀⢁⠀⠀⠈⠀⠁⠀⠀⠀⢀⠁⠄⠈⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠈⠀⠂⠠⠀⠀⠀⠀⠀⡀⠐⎥
⎢⠀⠀⠀⠈⠁⢀⠀⠑⠀⢀⠁⢀⠀⠈⠀⠘⠌⠀⢀⠀⠄⠀⠂⡄⠄⠁⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⣂⠀⠀⠀⎥
⎢⠀⢀⠁⠐⢄⠀⠠⢆⠀⠀⠀⠀⠀⠀⠀⢀⠠⡀⠀⠀⠠⠀⠀⠀⠐⠀⠀⡀⠀⢀⠀⠀⠀⠈⠀⡀⠀⠀⠘⠀⎥
⎢⠀⠐⡀⠀⠀⢀⠀⠀⢀⢔⠈⢀⠀⠀⣐⠀⠀⠀⢀⠀⠀⠀⠀⠀⠐⠀⠄⢠⠀⠀⠀⠀⠀⠀⠀⡀⠀⠈⠣⡀⎥
⎢⠀⠁⠄⠀⠁⢀⠀⠀⠂⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⡀⠀⠀⡐⠐⠀⡀⠀⠀⠂⠀⠀⠀⢀⠀⠀⠄⠀⎥
⎢⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠑⠀⠀⠀⠀⠀⠀⠀⡠⠀⡀⠀⠀⠄⠀⠀⠠⠀⠀⠠⠀⠀⠀⠀⡠⠄⠀⠄⎥
⎢⠀⠀⠀⢀⣀⠀⠀⢀⠐⠘⠀⠀⠀⠀⠕⢅⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀⢀⠐⡀⠀⠈⠀⠂⠀⠀⢀⠀⠃⠀⠄⎥
⎢⠀⠀⠁⠄⠂⠁⠀⠢⠀⠀⠀⠀⠀⠀⠀⠀⠛⢄⢸⠘⠀⠀⠀⠀⠄⠈⠁⠀⠀⠨⠀⠀⢀⠀⠀⠨⠀⠀⠈⠀⎥
⎢⠀⠄⠂⠀⠀⠐⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⣒⠒⠁⠀⢠⠀⠀⠀⠐⠀⠀⠀⢁⠀⠐⠈⠀⠀⠂⠀⠂⡀⠀⠀⎥
⎢⢀⠠⠀⠀⠀⠁⠀⠂⠀⠀⠀⠂⠀⠊⠀⠀⠀⠀⠀⠒⠑⠀⠀⠀⠀⠂⢀⠁⡂⠀⠀⢀⠀⠀⠀⠀⠁⠂⡊⠂⎥
⎢⢀⠀⢀⠀⠈⠤⠀⠀⠀⠀⠀⠈⠀⠈⠀⠠⠀⠀⠀⠀⠀⠀⠁⢀⠅⠀⢀⠀⠀⠀⢀⠀⡁⠀⠀⠀⠀⠠⠀⠀⎥
⎢⠠⠀⠀⠀⠄⠁⠐⠀⠐⠀⢀⠠⠀⠄⠀⠀⡀⠁⠐⠀⠠⠀⠁⠁⠁⠀⢀⠄⠀⠉⠀⠃⠀⠀⠀⠀⠠⢠⠂⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠠⠀⣁⠐⠀⠀⠀⢀⠐⠁⠀⠀⠀⠄⠐⠀⠐⠀⠔⡑⠌⠀⠀⠀⠀⠢⢉⠀⠀⠄⠀⠀⠄⎥
⎢⠀⠀⡀⠀⠂⠀⠀⢀⠀⠀⠀⠈⠀⠂⠀⠈⡀⡀⠁⠐⠈⠈⠀⠀⡄⠀⠀⠀⠄⠅⠀⡀⠀⠀⠀⠠⠀⠰⠀⠂⎥
⎢⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀⡀⠂⠀⠀⠀⡐⠀⠀⢀⠀⠐⠤⠀⠀⠀⠀⠠⠀⢀⠀⠀⠀⠀⠀⠀⠒⢀⎥
⎢⠈⠈⠀⠂⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠐⠀⠀⠀⠀⠁⠈⠀⠀⡌⢂⠀⠀⠀⠀⠀⠄⠈⠀⠀⡀⠀⠀⎥
⎢⠀⠔⠀⠀⠀⠀⠀⠠⠀⠠⠀⢀⠀⠀⠀⢀⡀⡀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀⠂⠀⠡⠂⠀⠀⡠⠀⎥
⎢⠁⡀⠀⠀⠈⠘⠀⠀⡀⠀⠀⠀⠀⠎⠤⠀⠀⠀⠈⠠⠡⠀⠀⡀⠀⣂⠀⠁⢀⡀⠀⠀⠀⠠⠀⠀⠰⠆⠌⠀⎥
⎣⠐⠀⢀⠈⠀⠀⠒⠀⠉⠢⠀⠁⠀⠄⠀⠄⠂⠀⠀⠀⠪⠈⠀⠀⠈⠀⠀⠄⠠⠀⠘⢀⠀⠀⠀⠊⠂⠁⠐⢀⎦

julia> A_shuffled[res.ordering, res.ordering] # A near-optimal reordering of `A_shuffled`
235×235 SparseMatrixCSC{Float64, Int64} with 445 stored entries:
⎡⠁⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠁⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠁⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠑⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢛⢔⢤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠳⣿⣿⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠚⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠡⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠫⣦⣤⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⢿⣷⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠯⣧⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠯⣧⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠺⢆⡄⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢻⣲⣄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢻⣶⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣢⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠱⣦⎦

julia> bandwidth(A)
9

julia> bandwidth(A_shuffled) # Much larger after shuffling
226

julia> res # Gets very close to the original bandwidth
Results of Bandwidth Minimization Algorithm
 * Algorithm: Reverse Cuthill–McKee
 * Approach: heuristic
 * Minimum Bandwidth: 9
 * Original Bandwidth: 226
 * Matrix Size: 235×235
```

# Notes
Note that the `node_selector` field must be of the form
`(A::AbstractMatrix{Bool}) -> Integer` (i.e., it must take in an boolean matrix and return
an integer). If this is not the case, an `ArgumentError` is thrown upon construction.

See also the documentation for [`CuthillMcKee`](@ref)—the original (non-reversed) algorithm.
(Indeed, the reverse Cuthill–McKee method of `_bool_minimal_band_ordering` is merely a
wrapper around the Cuthill–McKee method.)
"""
struct ReverseCuthillMcKee <: HeuristicSolver
    node_selector::Function

    function ReverseCuthillMcKee(node_selector::Function=DEFAULT_SELECTOR)
        _assert_valid_node_selector(node_selector)
        return new(node_selector)
    end
end

Base.summary(::ReverseCuthillMcKee) = "Reverse Cuthill–McKee"

_requires_symmetry(::ReverseCuthillMcKee) = true

#= We take advantage of the laziness of `Iterators.map` and `Iterators.flatmap` to avoid
allocating `component_orderings` or individual `component[component_ordering]` arrays.
(Indeed, the only allocations performed here are those performed by `_connected_components`
, individual `_cm_connected_ordering` calls, and `collect` at the very end.) =#
function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::CuthillMcKee)
    node_selector = solver.node_selector
    components = _connected_components(A)

    component_orderings = Iterators.map(
        component -> _cm_connected_ordering(view(A, component, component), node_selector),
        components,
    )

    return collect(
        Iterators.flatmap(
            ((component, component_ordering),) -> component[component_ordering],
            zip(components, component_orderings),
        ),
    )
end

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::ReverseCuthillMcKee)
    return reverse!(_bool_minimal_band_ordering(A, CuthillMcKee(solver.node_selector)))
end

# Cuthill–McKee searches each connected component independently
function _cm_connected_ordering(A::AbstractMatrix{Bool}, node_selector::Function)
    n = size(A, 1)

    ordering = Vector{Int}(undef, n)
    visited = falses(n)
    degrees = vec(sum(A; dims=1))
    queue = Queue{Int}()

    start = node_selector(A)
    visited[start] = true
    enqueue!(queue, start)

    for i in 1:n
        parent = dequeue!(queue)
        ordering[i] = parent

        unvisited = filter!(node -> !visited[node], findall(view(A, :, parent)))
        sort!(unvisited; by=node -> degrees[node])

        visited[unvisited] .= true
        foreach(neighbor -> enqueue!(queue, neighbor), unvisited)
    end

    return ordering
end
