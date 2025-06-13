# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    ReverseCuthillMcKee <: HeuristicSolver <: AbstractSolver

The *reverse Cuthill–McKee algorithm* is a variant of the *Cuthill–McKee algorithm*—a
heuristic method for minimizing the bandwidth of a symmetric matrix ``A``. Cuthill–McKee
considers the graph ``G(A)`` whose adjacency matrix is ``A`` (ignoring self-loops) and
performs a breadth-first search of each connected component of ``G(A)``, starting from a
low-degree node then visiting its neighbors in order of increasing degree. Particularly
effective when ``A`` is sparse, this heuristic typically produces an ordering which induces
a matrix bandwidth either equal to or very close to the true minimum
[CM69; pp. 157--58](@cite). The reverse Cuthill–McKee algorithm simply reverses the ordering
produced by application of Cuthill–McKee; it was found in [Geo71; pp. 114--15](@cite) that
this tends to induce an even more optimal bandwidth.

We also extend the algorithm to work more generally when ``A`` is not symmetric by applying
it to ``A + Aᵀ`` instead, as suggested in [RS06; p. 808](@cite). This approach still tends
to produce a fairly good ordering, but it is not guaranteed to be as optimal as directly
applying reverse Cuthill–McKee to a symmetric input.

# Performance
Given an ``n×n`` input matrix ``A``, our implementation of Cuthill–McKee runs in ``O(n^2)``
time, where ``n`` is the number of rows/columns of the input matrix.

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

# Examples
Reverse Cuthill–McKee finds an optimal ordering for an asymmetric ``45×45`` matrix whose
rows and columns have been shuffled:
```jldoctest
julia> using Random

julia> Random.seed!(87);

julia> (n, k) = (45, 4);

julia> A = random_banded_matrix(n, k);

julia> perm = randperm(n);

julia> A_shuffled = A[perm, perm];

julia> iszero.(A) != iszero.(A') # Proof that the algorithm works for asymmetric input
true

julia> bandwidth(A)
4

julia> bandwidth(A_shuffled) # Much larger after shuffling
44

julia> res = minimize_bandwidth(A_shuffled, ReverseCuthillMcKee()) # Finds the true minimum
Results of Matrix Bandwidth Minimization
 * Algorithm: Reverse Cuthill–McKee algorithm
 * Approach: Heuristic
 * Minimum bandwidth: 4
 * Original bandwidth: 44
 * Matrix size: 45×45
```

Reverse Cuthill–McKee finds a near-optimal ordering for an asymmetric ``251×251`` matrix
with multiple (separate) connected components whose rows and columns have been shuffled:
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

julia> res = minimize_bandwidth(A_shuffled, ReverseCuthillMcKee());

julia> A # The original matrix
251×251 SparseMatrixCSC{Float64, Int64} with 617 stored entries:
⎡⢛⣷⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠈⢿⣇⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠘⠻⣦⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠉⢿⣶⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠘⠹⣢⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢿⣷⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠬⣣⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠿⣵⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢛⣧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢾⣳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢻⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⡓⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⡥⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠾⢇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⡆⡀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠿⣡⡄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠻⣤⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠳⣆⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣆⎦

julia> A_shuffled # A far-from-optimal ordering of `A`
251×251 SparseMatrixCSC{Float64, Int64} with 617 stored entries:
⎡⠑⢅⠀⠀⡠⠁⡀⠀⠨⢀⠀⠀⠄⠀⠀⡂⠀⡁⠒⢄⠀⠂⢀⠀⠀⠐⣁⢀⠀⢂⠠⠁⡠⠢⢀⡂⠀⠀⢀⡀⎤
⎢⠀⠀⠩⢂⠀⠀⠀⣀⠈⠈⠁⠀⠀⠠⠈⠀⠀⠀⠄⠀⠀⠒⠀⠀⠀⠀⠀⠈⠀⠁⠈⡀⠀⠀⠠⡀⠀⠀⠄⠀⎥
⎢⡀⠦⠀⠀⠑⠄⠀⠀⠠⠀⠠⠄⠀⠀⢀⠌⠀⠀⠀⠈⢠⠀⠁⡀⠢⠀⠀⠠⡁⠂⡔⡀⠁⠀⢂⠀⠠⢁⠀⠁⎥
⎢⠠⠈⠄⣠⠐⠀⠐⢄⡁⠀⠀⡀⢰⠀⠁⠀⠈⠀⠀⠀⢀⠂⡀⠀⢀⠀⠁⠂⢂⠈⠀⢀⠀⠀⠀⠈⠀⠈⡀⡠⎥
⎢⠀⠀⠀⠄⠀⡈⠀⠀⠵⠂⠀⠀⠀⢠⠀⠂⠒⠀⠀⠀⠐⠠⠀⠀⢄⠀⣐⠆⠀⠐⠀⠠⠀⠀⠨⠄⠦⡀⠀⡀⎥
⎢⠠⠀⢀⠀⠀⠀⠀⠠⠀⠈⡐⢌⠀⠂⠂⠈⠀⠔⠀⠁⣀⠀⢀⠀⠠⠀⠀⠐⠀⠀⠐⠁⠀⠄⠀⠂⠀⠀⠂⢀⎥
⎢⢨⢁⠀⠠⠀⠀⠀⠄⡀⡀⠀⠁⠂⠂⠁⠐⢀⠀⠀⢁⠐⠘⢂⠀⠀⠀⡀⠀⠔⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⎥
⎢⠢⠂⠂⡀⠀⠔⠀⢂⠀⠀⢨⠀⢀⠀⠐⠐⠀⠀⠀⠂⠐⠄⠈⠀⡀⠂⠁⠤⢀⠂⢀⠀⡀⢀⠀⠀⠀⠀⠀⠀⎥
⎢⠂⠀⠀⠀⠀⠀⠂⠀⠄⠂⢀⠁⠀⠰⡀⠀⢕⢐⠀⠀⠀⠀⠂⡀⠠⠄⠀⠁⠂⠈⠀⠄⠚⠀⠐⠀⠀⠀⢀⠁⎥
⎢⠀⠄⠀⠁⠁⠀⡀⠀⠀⠀⠁⠀⠄⠘⡀⠠⡊⠀⠒⠄⠐⠂⠀⡂⠁⠀⠂⠀⠀⠀⠀⢄⠁⠀⠀⠀⠀⠀⠠⡀⎥
⎢⠈⠘⠀⠂⠀⠀⠀⠐⠀⠃⠀⠀⠀⡀⠀⠄⠁⠀⠄⠀⠎⢑⠀⢂⠀⠈⠁⠀⠀⢀⠀⢀⠀⢀⠰⡁⠀⠢⠀⠀⎥
⎢⢄⠀⡀⢀⠃⠀⠁⡈⠁⠀⠠⠀⠀⠄⠄⠠⠠⠀⠄⠀⠒⠁⠱⠆⠀⠀⠀⠀⡀⠀⠀⠄⠁⠁⢀⠀⠀⠊⠁⠀⎥
⎢⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠐⠂⠀⠀⡄⠀⠠⠀⠂⠀⠀⠀⠘⠈⠔⠀⠀⠀⠀⠀⠀⠁⠐⢁⢀⠈⡀⠠⠐⎥
⎢⠀⠐⡀⠀⠠⢠⠠⠀⠒⢀⠀⠀⠀⠈⠀⡄⠀⠄⠀⠀⠈⠀⠁⠀⠀⠀⢶⢐⠄⠀⠀⠤⠀⠀⠐⡐⡘⠀⡁⠄⎥
⎢⠁⠀⠒⠀⠀⡀⠐⠁⡂⠀⠀⠀⠔⢀⠈⠀⠈⠀⠀⠒⠀⠄⠂⠀⠒⠀⠀⠀⠀⠄⡀⢐⠐⢄⡄⠀⢀⠀⡄⠈⎥
⎢⠀⠀⠀⠀⠑⠀⠀⠀⠀⠀⠄⠀⠀⠀⠐⠸⠀⠄⠀⢄⠀⠀⡀⡄⠁⠆⢀⠀⠠⠀⠑⠀⠅⢁⠀⠀⠨⠄⠀⠐⎥
⎢⠠⡂⠀⠀⠀⠂⡀⡀⠀⠁⠀⠀⠀⡠⠀⠀⠐⠀⠀⡉⠀⠀⠀⢀⠐⠀⠉⠀⠀⢄⠀⠁⠐⢔⢀⡀⢂⠀⠀⠀⎥
⎢⢀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠐⠀⠐⠀⠀⢄⠠⠀⠀⠀⠰⢀⡐⠂⢐⠀⠁⠀⠁⠀⠀⢀⠐⠅⢅⠄⢀⠂⠈⎥
⎢⠀⠀⠀⠀⠢⢀⠠⠀⠈⠂⠀⠀⡀⠐⡀⠀⠑⠐⠂⠈⠪⠐⠈⠂⠀⠀⠂⠈⠀⠐⠀⠄⠀⠄⠀⠐⠑⢀⠐⠀⎥
⎣⠀⠀⡀⠁⠐⠀⠀⠈⠀⠉⠉⢀⡀⠀⠁⠀⡈⠀⠀⠂⠀⡀⠀⠀⠀⠂⠀⡀⡄⠀⠁⠂⠁⠀⠀⠀⠀⠀⠑⢆⎦

julia> A_shuffled[res.ordering, res.ordering] # A near-optimal reordering of `A_shuffled`
251×251 SparseMatrixCSC{Float64, Int64} with 617 stored entries:
⎡⠑⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⢿⣷⣠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠈⠹⢷⣷⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠈⢻⡶⣆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⣾⣧⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣧⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢛⣟⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠓⢿⣲⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢿⣷⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠻⡦⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠿⣧⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢿⣶⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠫⣶⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠿⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⢻⣴⡆⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢱⢶⡄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⣦⣄⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠾⣷⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠯⣆⎦

julia> iszero.(A) != iszero.(A') # Proof that the algorithm works for asymmetric input
true

julia> bandwidth(A)
7

julia> bandwidth(A_shuffled) # Much larger after shuffling
239

julia> res # Gets very close to the true minimum
Results of Matrix Bandwidth Minimization
 * Algorithm: Reverse Cuthill–McKee algorithm
 * Approach: Heuristic
 * Minimum bandwidth: 10
 * Original bandwidth: 239
 * Matrix size: 251×251
```

# Notes
Note that the `node_selector` field must be of the form
`(A::AbstractMatrix{Bool}) -> Integer` (i.e., it must take in an boolean matrix and return
an integer). If this is not the case, an `ArgumentError` is thrown upon construction.

See also the documentation for supertypes [`HeuristicSolver`](@ref) and
[`AbstractSolver`](@ref), as well as [`CuthillMcKee`](@ref) for the original non-reversed
algorithm. (Indeed, the reverse Cuthill–McKee method of `_bool_minimal_band_ordering` is
merely a wrapper around the Cuthill–McKee method.)
"""
struct ReverseCuthillMcKee <: HeuristicSolver
    node_selector::Function

    function ReverseCuthillMcKee(node_selector::Function=DEFAULT_SELECTOR)
        _assert_valid_node_selector(node_selector)
        return new(node_selector)
    end
end

Base.summary(::ReverseCuthillMcKee) = "Reverse Cuthill–McKee algorithm"

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::ReverseCuthillMcKee)
    return reverse!(_bool_minimal_band_ordering(A, CuthillMcKee(solver.node_selector)))
end
