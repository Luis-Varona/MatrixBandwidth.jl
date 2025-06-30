# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CuthillMcKee <: HeuristicSolver <: AbstractSolver

The *Cuthill–McKee algorithm* is a heuristic method for minimizing the bandwidth of a
symmetric matrix ``A``. It considers the graph ``G(A)`` whose adjacency matrix is ``A``
(ignoring self-loops) and performs a breadth-first search of each connected component of
``G(A)``, starting from a low-degree node then visiting its neighbors in order of increasing
degree. Particularly effective when ``A`` is sparse, this heuristic typically produces an
ordering which induces a matrix bandwidth either equal to or very close to the true minimum
[CM69; pp. 157--58](@cite).

We also extend the algorithm to work more generally when ``A`` is not symmetric by applying
it to ``A + Aᵀ`` instead, as suggested in [RS06; p. 808](@cite). This approach still tends
to produce a fairly good ordering, but it is not guaranteed to be as optimal as directly
applying Cuthill–McKee to a symmetric input.

# Fields
- `node_selector::Function`: a function that selects a node from some connected component of
    the input matrix from which to start the breadth-first search. If no custom heuristic is
    specified, this field defaults to [`pseudo_peripheral_node`](@ref), which picks a node
    "farthest" from the others in the component (not necessarily the lowest-degree node).

# Performance
Given an ``n×n`` input matrix ``A``, our implementation of Cuthill–McKee runs in ``O(n^2)``
time, where ``n`` is the number of rows/columns of the input matrix.

[CG80](@cite) provide a linear-time implementation in the number of nonzero entries of
``A``, which is still quadratic when ``A`` is dense but often much faster when dealing with
sparse matrices. However, this would require that ``A`` be stored as a graph or a sparse
matrix, which runs counter to our desire to provide a bandwidth minimization API for all
`AbstractMatrix{<:Number}` types, including dense matrices. (In the future, however, we may
indeed consider supporting this more performant implementation for sparse matrices.)

It was found in [Geo71; pp. 114--15](@cite) that reversing the ordering produced by
Cuthill–McKee tends to induce a more optimal *matrix profile* (a measure of how far, on
average, nonzero entries are from the diagonal). This so-called *reverse Cuthill–McKee*
variant is preferred in almost all cases—see [`ReverseCuthillMcKee`](@ref) and the
associated method of `_bool_minimal_band_ordering` for our implementation.

# Examples
Cuthill–McKee finds an optimal ordering for an asymmetric ``35×35`` matrix whose rows and
columns have been shuffled:
```jldoctest
julia> using Random

julia> Random.seed!(13);

julia> (n, k) = (35, 3);

julia> A = random_banded_matrix(n, k);

julia> perm = randperm(n);

julia> A_shuffled = A[perm, perm];

julia> iszero.(A) != iszero.(A') # Proof that the algorithm works for asymmetric input
true

julia> bandwidth(A)
3

julia> bandwidth(A_shuffled) # Much larger after shuffling
31

julia> res = minimize_bandwidth(A_shuffled, CuthillMcKee()) # Finds the true minimum
Results of Matrix Bandwidth Minimization
 * Algorithm: Cuthill–McKee algorithm
 * Approach: Heuristic
 * Minimum bandwidth: 3
 * Original bandwidth: 31
 * Matrix size: 35×35
```

Cuthill–McKee finds a near-optimal ordering for an asymmetric ``183×183`` matrix with
multiple (separate) connected components whose rows and columns have been shuffled:
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

julia> res = minimize_bandwidth(A_shuffled, CuthillMcKee());

julia> A # The original matrix
183×183 SparseMatrixCSC{Float64, Int64} with 408 stored entries:
⎡⣜⣹⡤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠪⣿⣭⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠁⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠑⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢱⢆⢠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠳⢽⡇⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⠽⡺⠦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠓⣷⣇⣂⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠣⣏⣾⣂⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠲⣻⣾⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠋⡏⢷⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⠾⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⠄⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⢽⣟⣀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⠾⠶⡤⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠊⠾⡻⣦⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠃⡟⣵⢄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠱⢹⣣⎦

julia> A_shuffled # A far-from-optimal ordering of `A`
183×183 SparseMatrixCSC{Float64, Int64} with 408 stored entries:
⎡⠁⢄⡂⠀⠀⢀⠀⠀⠀⠂⠂⢀⠀⢀⠀⠀⠀⠀⠀⢀⠐⠠⠂⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⎤
⎢⠀⠈⠁⠤⠀⠄⠀⠀⠀⠀⠠⠀⢀⠀⠀⠁⠡⠀⠀⠀⠈⠀⠀⡀⠀⡀⠒⠀⠀⠘⠀⠀⠀⠀⡀⠀⠀⠀⠀⡄⎥
⎢⠀⢠⠀⠀⠁⢄⠈⠀⠀⠀⠀⠀⠄⠀⠀⠀⠠⠀⠐⠁⠀⠠⠁⠠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠤⠀⠀⎥
⎢⠀⠀⠀⠀⠢⠀⠁⠅⠀⠀⠐⢐⠀⠂⠀⠀⠀⠀⠐⠐⠀⠐⠀⠀⠀⠀⠀⠀⠀⠄⠠⠀⠀⠐⠐⠀⠀⠀⠀⠂⎥
⎢⠀⠀⠀⠀⠉⢀⢀⠁⠁⢀⠀⡈⠁⠐⠀⠀⠐⠀⠈⠀⠀⡀⠀⠀⠀⠀⠈⠀⠀⠈⠀⠀⠂⠀⠀⠐⠐⠀⠐⠀⎥
⎢⠀⢀⢀⠀⡀⠁⠀⠀⠈⠀⠑⣀⡀⠀⠀⠀⠐⠀⠀⠈⠀⠀⠈⠇⠀⠀⠐⠀⠀⠀⠀⠀⠠⠀⠀⠄⠁⠂⠁⠐⎥
⎢⠀⠀⠀⠐⠀⠈⠠⠀⠀⠈⡀⠀⢐⠌⠀⡁⡀⠐⠀⠰⠀⢀⠀⠄⠀⠄⡀⠀⠀⠀⠐⠔⡁⠀⠀⠀⠀⠀⡀⠄⎥
⎢⠀⠀⠀⠠⠄⠀⠀⠀⠀⠀⠁⠀⠀⠀⠐⢈⠤⡀⠀⠄⠈⠀⠀⠀⠀⠄⠀⠁⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⢀⠀⎥
⎢⠈⠀⠁⢀⠀⠀⠀⠀⠐⠀⢀⠀⠀⠂⠀⠠⠑⠀⠀⠀⠈⠂⠀⠂⠐⡀⢀⠤⠂⠀⠐⠀⠀⠀⠐⡠⠀⠂⠀⠀⎥
⎢⢀⠀⠀⢀⠆⡠⠂⠀⠂⠀⠀⠀⠄⠀⠀⠂⠀⠁⠀⢀⠀⠔⠀⠀⠀⡀⠀⠀⠂⢂⠀⠀⠀⢀⠀⠀⠀⢀⠀⠃⎥
⎢⠐⡀⠂⡀⠀⠀⢀⠀⢄⠠⡀⠀⠀⠀⠀⠀⠀⠈⡁⡀⠀⢤⠒⡀⠀⠀⠀⠈⠀⠀⠠⠀⠀⠐⠀⠁⠀⠀⡀⠀⎥
⎢⠀⠀⠀⠠⠠⠀⠀⠀⠁⠀⠂⠈⠀⠀⠀⡀⠁⢤⠄⠈⠈⠠⠀⠄⠈⠀⠀⠠⠀⠀⠀⠀⠈⠀⠠⠌⡀⠈⠠⠀⎥
⎢⠀⠐⠀⠀⠀⠀⠠⠀⠀⠀⠠⠀⠀⠀⠀⠀⠄⠁⠀⠲⠀⠀⠀⠀⠐⠀⠂⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠠⠀⠀⠐⠀⠀⠀⠈⠀⠀⠀⠀⠀⡀⠀⠀⡀⠀⠀⡄⠀⠀⠐⠄⠀⠀⠐⠀⠊⠀⠐⡢⠈⠀⠀⢀⎥
⎢⠀⠀⡠⠀⠁⡀⠌⠄⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠌⠀⠀⠀⠀⠀⠀⠀⣀⠀⠑⠀⢀⡀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⡐⢂⠀⠅⠀⠂⠀⠀⠀⠀⠀⠀⠄⠀⠀⢀⠀⠠⡀⠂⠄⡀⠀⠀⠂⠀⠀⠘⠑⠀⠀⠀⠀⠀⠠⠢⠀⠀⎥
⎢⠀⠀⠂⠀⠠⢀⠀⠀⠀⠀⠀⢀⠀⠄⠀⠄⠀⠄⠄⠠⢀⠀⠀⠀⠀⠀⠁⠠⠀⢘⠀⠠⠐⠀⠁⠀⠀⢀⠠⠀⎥
⎢⠀⢂⠀⠀⠀⠀⠀⠀⠐⠀⠀⡀⡀⠀⠀⢀⠀⠈⠀⠀⠀⠀⡀⠂⠐⠀⠸⠨⡀⠀⠀⠀⢀⡀⠱⠄⠈⠀⠀⠈⎥
⎢⢀⠐⠀⠠⠀⠂⠈⠀⠀⠀⠀⠀⠂⠡⠀⠄⠀⠀⠀⢠⠀⠠⡀⠀⠀⠀⠀⠀⠀⢐⠀⠀⠀⠐⠂⠀⠁⢐⠂⠀⎥
⎣⠀⠀⠀⠄⠀⢀⠀⠀⠀⠀⠁⠀⠠⠀⠀⠤⠀⠂⠤⠀⠀⠀⠀⠂⠈⠀⠀⠐⠀⠀⠀⠀⠀⠀⠠⠠⠈⠄⠀⠄⎦

julia> A_shuffled[res.ordering, res.ordering] # A near-optimal reordering of `A_shuffled`
183×183 SparseMatrixCSC{Float64, Int64} with 408 stored entries:
⎡⢱⣶⣤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠹⡵⣉⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠢⡤⡤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⢳⡺⣺⠆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠛⡗⡻⣤⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⡎⣽⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⣟⣎⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠘⠬⢅⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢳⣓⠆⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢼⣥⣒⢂⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢲⣫⣾⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠳⣜⡼⡤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠃⣽⠳⣲⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠐⠛⢄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⣿⣘⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠁⠃⢦⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢄⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠓⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠄⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠀⎦

julia> iszero.(A) != iszero.(A') # Proof that the algorithm works for asymmetric input
true

julia> bandwidth(A)
7

julia> bandwidth(A_shuffled) # Much larger after shuffling
170

julia> res # Gets very close to the true minimum
Results of Matrix Bandwidth Minimization
 * Algorithm: Cuthill–McKee algorithm
 * Approach: Heuristic
 * Minimum bandwidth: 10
 * Original bandwidth: 170
 * Matrix size: 183×183
```

# Notes
Note that the `node_selector` field must be of the form
`(A::AbstractMatrix{Bool}) -> Integer` (i.e., it must take in an boolean matrix and return
an integer). If this is not the case, an `ArgumentError` is thrown upon construction.

See also the documentation for supertypes [`HeuristicSolver`](@ref) and
[`AbstractSolver`](@ref).
"""
struct CuthillMcKee <: HeuristicSolver
    node_selector::Function

    function CuthillMcKee(node_selector::Function=DEFAULT_SELECTOR)
        _assert_valid_node_selector(node_selector)
        return new(node_selector)
    end
end

Base.summary(::CuthillMcKee) = "Cuthill–McKee algorithm"

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::CuthillMcKee)
    n = size(A, 1)
    A_sym = _symmetrize(A)
    A_sym[1:(n + 1):end] .= false

    node_selector = solver.node_selector
    components = _connected_components(A_sym)
    ordering = Vector{Int}(undef, n)
    k = 1

    for component in components
        submatrix = A_sym[component, component]
        component_ordering = _connected_cuthill_mckee_ordering(submatrix, node_selector)

        component_size = length(component)
        ordering[k:(k + component_size - 1)] = component[component_ordering]
        k += component_size
    end

    return ordering
end

# Cuthill–McKee performs a breadth-first search on each connected component independently
function _connected_cuthill_mckee_ordering(A::AbstractMatrix{Bool}, node_selector::Function)
    n = size(A, 1)
    ordering = Vector{Int}(undef, n)

    start = node_selector(A)
    degrees = vec(sum(A; dims=1))
    visited = Set(start)
    queue = Queue{Int}()
    enqueue!(queue, start)

    for i in 1:n
        parent = dequeue!(queue)
        ordering[i] = parent

        neighbors = findall(A[:, parent])
        unvisited = filter(!in(visited), neighbors)
        sort!(unvisited; by=node -> degrees[node])

        union!(visited, unvisited)
        foreach(neighbor -> enqueue!(queue, neighbor), unvisited)
    end

    return ordering
end

# Find the indices of all connected components in an adjacency matrix
function _connected_components(A::AbstractMatrix{Bool})
    n = size(A, 1)
    visited = falses(n)
    queue = Queue{Int}()
    components = Vector{Int}[]

    for i in 1:n
        if !visited[i]
            visited[i] = true
            enqueue!(queue, i)
            component = Int[]

            while !isempty(queue)
                u = dequeue!(queue)
                push!(component, u)

                for v in 1:n
                    if A[u, v] && !visited[v]
                        visited[v] = true
                        enqueue!(queue, v)
                    end
                end
            end

            push!(components, component)
        end
    end

    return components
end
