# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    GibbsPooleStockmeyer <: HeuristicSolver <: AbstractSolver <: AbstractAlgorithm

The *Gibbs–Poole–Stockmeyer algorithm* is a heuristic method for minimizing the bandwidth of
a structurally symmetric matrix ``A``. It considers the graph ``G(A)`` whose adjacency
matrix is ``A`` (ignoring weights and self-loops) and builds an ordering by identifying a
pair of "endpoints" in the graph far from each other, constructing sets of levels from
these endpoints, and merging these level structures in such a way that minimizes the size of
the largest level in the final combined structure. Based on the classical reverse
Cuthill–McKee algorithm [Geo71], this heuristic typically produces an ordering which induces
a matrix bandwidth either equal to or very close to the true minimum, with improvements in
bandwidths over reverse Cuthill–McKee more noticeable once input size exceeds ``400×400`` or
so [GPS76, pp. 246--47].

Whereas the original paper outlined a strategy for conditionally reversing the orderings of
individual "connected components" (treating the input matrix ``A`` as an undirected graph)
[GPS76, p. 241], this implementation instead reverses the entire final ordering in every
case, similarly to [`ReverseCuthillMcKee`](@ref). Conditional reversals are not only more
complex to implement but also slightly more time-consuming, with the only benefit being a
marginally smaller *matrix profile* (a measure of how far, on average, nonzero entries are
from the diagonal; see also [`MatrixBandwidth.profile`](@ref)). Since such reversal
strategies do not affect matrix bandwidth (the primary focus of this package), we opt
instead for the simpler unconditional reversal.

As noted above, the Gibbs–Poole–Stockmeyer algorithm requires structurally symmetric input
(that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for
``1 ≤ i, j ≤ n``).

# Fields
- `node_finder::Function`: a function that selects a node from some connected component of
    the input matrix from which to start the breadth-first search. If no custom heuristic is
    specified, this field defaults to [`bi_criteria_node_finder`](@ref), which picks a node
    "farthest" from the others in the component (not necessarily the lowest-degree node).

# Supertype Hierarchy
`GibbsPooleStockmeyer` <: [`HeuristicSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix, the Gibbs–Poole–Stockmeyer algorithm runs in ``O(n²)`` time.

[Lew82] provides a notably faster and more memory-efficient implementation, relying on
sparse storage of the input matrix. However, this would run counter to our desire to provide
a bandwidth minimization API for all `AbstractMatrix{<:Number}` types, including dense
matrices. (In the future, however, we may indeed consider supporting this more performant
implementation for sparse matrices.)

On that note, Gibbs–Poole–Stockmeyer has been found to take considerably less time than
reverse Cuthill–McKee when matrices are stored in sparse format [GPS76, pp. 246--47]. The
dense-matrix implementations of both algorithms in this package, however, result in reverse
Cuthill–McKee consistently outperforming Gibbs–Poole–Stockmeyer in terms of runtime
(although Gibbs–Poole–Stockmeyer still typically produces lower-bandwidth orderings for
larger matrices). This further motivates the desire to implement a sparse version of both
algorithms in the future.

# Examples
In the following examples, [`MatrixBandwidth.random_banded_matrix`](@ref) is used to
generate random matrices with minimum bandwidth *close to* ``k``. In some cases, however,
the true minimum bandwidth up to symmetric permutation may be even less than ``k``, making
it hard to verify whether Gibbs–Poole–Stockmeyer finds a truly optimal ordering or simply a
near-optimal one. Nevertheless, the results are still very good in practice.

Gibbs–Poole–Stockmeyer finds a good ordering for a ``40×40`` matrix:
```jldoctest
julia> using Random

julia> Random.seed!(561);

julia> (n, k) = (40, 7);

julia> A = MatrixBandwidth.random_banded_matrix(n, k);

julia> perm = randperm(n);

julia> A_shuffled = A[perm, perm];

julia> bandwidth(A)
7

julia> bandwidth(A_shuffled)
37

julia> minimize_bandwidth(A_shuffled, Minimization.GibbsPooleStockmeyer())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Gibbs–Poole–Stockmeyer
 * Approach: heuristic
 * Minimum Bandwidth: 7
 * Original Bandwidth: 37
 * Matrix Size: 40×40
```

Gibbs–Poole–Stockmeyer finds a good ordering for a ``748×748`` matrix with multiple
(separate) connected components:
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(271828);

julia> (max_cc_size, max_band, p, num_ccs) = (120, 13, 0.3, 11);

julia> components = Vector{SparseMatrixCSC{Float64, Int64}}(undef, num_ccs);

julia> for i in 1:num_ccs # Some components may themselves be disconnected
           cc_size = rand(0:max_cc_size);
           cc_band = rand(1:min(max_band, cc_size - 1));
           components[i] = sparse(
               MatrixBandwidth.random_banded_matrix(cc_size, cc_band; p=p)
           );
       end

julia> A = blockdiag(components...); # `A` has least 8 connected components

julia> perm = randperm(sum(map(cc -> size(cc, 1), components)));

julia> A_shuffled = A[perm, perm];

julia> res = minimize_bandwidth(A_shuffled, Minimization.GibbsPooleStockmeyer());

julia> A # The original matrix
748×748 SparseMatrixCSC{Float64, Int64} with 2526 stored entries:
⎡⢿⣷⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠘⠿⣧⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠘⢿⣷⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠉⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⣤⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⎦

julia> A_shuffled # A far-from-optimal ordering of `A`
748×748 SparseMatrixCSC{Float64, Int64} with 2526 stored entries:
⎡⠰⣦⢎⢪⢐⠆⣗⣔⠆⠀⠀⠠⢃⡦⠵⠸⡐⢌⠴⠤⣤⢅⢒⠰⡢⡄⡰⣄⢂⢊⠎⠀⢝⡀⣼⠤⠅⢒⠰⠢⎤
⎢⡪⣑⡋⢌⠈⢀⣖⣠⢠⡀⡂⠀⠐⠤⣠⠰⣃⠀⢊⢃⡨⡇⡝⠢⣀⠈⢬⢁⠽⡼⠐⡦⠤⠐⠤⠚⠪⢠⠦⢈⎥
⎢⠰⠔⠂⢀⡑⢌⠐⠀⠐⡒⠤⡅⠂⡱⡕⠈⡑⠹⠱⣨⢓⣑⠀⠂⠀⢬⢂⣅⡤⢑⠢⢁⠀⢒⡑⡬⡆⠚⠀⡀⎥
⎢⢙⢽⠘⣹⠐⠀⢱⣶⠂⠂⠈⡵⢰⡤⠖⠰⡂⣠⢨⢰⣛⠀⠰⡌⡐⠄⠠⡦⠜⠃⠐⣦⢂⠄⠷⢅⣉⠰⠿⢴⎥
⎢⠈⠁⠀⠲⢰⠠⠨⠀⣑⣼⠩⡃⠂⢰⢁⠓⢐⣘⣂⢦⠂⡐⠔⢄⡨⣃⠦⢁⢈⠀⡂⠀⡅⠈⠀⢌⠀⡀⠵⠫⎥
⎢⠀⡀⠈⠈⠄⠧⢆⡤⠧⠢⠕⢅⡍⢔⡴⠀⡀⠆⡨⠈⣄⠲⠌⠳⢁⠐⠠⠴⠀⠄⠔⠺⠁⡉⢂⣭⠠⡠⠀⣉⎥
⎢⠩⡴⠐⡄⢌⡠⠐⡶⢈⣀⢃⢍⠋⢄⡖⢩⢌⣖⠈⢀⣴⣙⡀⢓⠁⠠⢬⢈⠅⡤⠅⢰⣁⠌⣌⠆⠸⢀⠒⣁⎥
⎢⣑⡃⢀⡚⡑⠉⢘⡁⢥⠐⠐⠋⡜⣉⠑⢄⠍⣢⣓⢉⠊⢖⠀⢐⡀⠲⡈⢑⠀⠊⠀⠛⠀⢃⠘⠀⢁⣒⢀⢁⎥
⎢⡐⢌⠉⠘⣕⡈⠈⣨⣐⢰⠠⠌⢢⢵⠣⣡⢕⣱⣂⢭⢎⠜⠀⠉⡂⣃⢐⠅⠒⠔⡂⠨⡂⢳⠍⢤⠰⡠⣩⡏⎥
⎢⠐⡇⠮⢐⡑⣢⢂⣒⠨⣜⡂⠊⠂⢀⡝⢘⡌⣜⣱⢞⠂⠾⠊⢔⠙⣢⢭⣹⠑⡑⠀⡉⡄⠁⠀⠴⣇⣃⠃⡑⎥
⎢⠄⢟⠦⠮⢝⢰⠛⠘⢈⠠⢠⡙⣔⢻⢪⢄⣊⠕⣨⡄⢕⣵⢐⠊⠊⣱⡠⣢⠀⠣⠀⣍⠔⠀⢝⢄⣌⡄⡌⠆⎥
⎢⢘⡐⠳⡉⠠⠀⡐⠦⠐⢅⢦⡁⢤⢈⢀⢀⡄⠀⢊⢄⡰⠐⠑⢄⠴⡅⠁⢆⡠⣄⠤⢐⠈⡀⠺⠂⠀⢀⡴⡀⎥
⎢⠈⠮⡀⠘⡀⣄⠐⠌⠦⢪⢁⠐⠁⡀⢠⡈⠬⢨⠳⣠⢎⣠⠔⠧⠛⢄⡀⡹⠡⠌⠃⢀⠣⠁⠉⢀⠑⠁⢦⡉⎥
⎢⠐⢮⠆⢓⠌⢴⠠⡦⠌⢃⢀⡆⡂⢓⢆⢈⠔⠔⣇⣳⠠⣪⠡⢄⣄⡨⠛⢄⠑⢐⠐⡠⠪⠃⢤⢁⡝⠐⡀⢎⎥
⎢⡨⢐⣓⡧⢄⢋⠶⠁⠂⠐⠀⠄⠁⡥⡠⠀⢘⠄⢕⠠⠤⡀⠀⢮⡁⠆⢑⢀⠕⣥⣑⠼⡀⡅⠝⠢⠠⠄⠀⠔⎥
⎢⠊⠁⠰⡤⠌⢂⠰⣤⠈⠈⣰⡁⢁⣁⣤⠀⡈⡈⡄⠠⡄⢤⢀⢃⠉⢀⠐⡠⣑⡜⠑⢄⠀⡼⠠⠲⢈⠠⢃⠁⎥
⎢⠓⠱⢀⠃⢠⢀⠈⠔⡁⠉⡅⠠⡁⠜⠤⢀⢬⣈⠄⠉⠐⠁⠂⠠⠍⠂⠮⠂⠄⠬⣀⡤⠵⢇⠈⠀⠣⠡⠩⡎⎥
⎢⠒⡟⣠⠃⡑⡬⠝⢇⡀⢄⡌⣴⠢⠝⠒⠀⠃⣅⢀⡄⠓⢕⠺⠂⠃⢀⠄⢓⠳⡁⢠⡂⠂⠀⠛⢄⡀⢒⠒⠂⎥
⎢⢡⢁⠊⣂⣨⠉⢃⡘⠀⠠⠀⡢⠒⢂⢡⢰⠐⡢⠭⢹⠂⠽⠀⢀⠕⠀⢓⠉⠀⠆⠂⡐⠍⡂⢠⢈⠁⣤⢀⡈⎥
⎣⠰⡂⡈⢃⠀⠠⢛⣇⡵⡃⡄⢠⠜⢠⠄⢐⡧⠾⢍⠠⠢⠍⠐⠫⡌⠳⡠⢌⢀⠄⠍⠐⡣⠦⠸⠀⡀⠰⡛⢌⎦

julia> A_shuffled[res.ordering, res.ordering] # A near-optimal reordering of `A_shuffled`
748×748 SparseMatrixCSC{Float64, Int64} with 2526 stored entries:
⎡⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠑⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠱⢆⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠱⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢻⣶⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣷⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⣤⣀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢻⣶⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⣤⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⎦

julia> bandwidth(A)
12

julia> bandwidth(A_shuffled) # Much larger after shuffling
731

julia> res # Gets very close to the original bandwidth
Results of Bandwidth Minimization Algorithm
 * Algorithm: Gibbs–Poole–Stockmeyer
 * Approach: heuristic
 * Minimum Bandwidth: 18
 * Original Bandwidth: 731
 * Matrix Size: 748×748
```

# Notes
Note that the `node_finder` field must be of the form `(A::AbstractMatrix{Bool}) -> Integer`
(i.e., it must take in an boolean matrix and return an integer). If this is not the case, an
`ArgumentError` is thrown upon construction.

# References
- [Geo71](@cite): J. A. George. *Computer Implementation of the Finite Element Method*.
    Ph.D. Thesis, Department of Computer Science, Stanford University (1971).
    https://apps.dtic.mil/sti/tr/pdf/AD0726171.pdf.
- [GPS76](@cite): N. E. Gibbs, W. G. Poole Jr. and P. K. Stockmeyer. *An Algorithm for
    Reducing the Bandwidth and Profile of a Sparse Matrix*. SIAM Journal on Numerical
    Analysis **13**, 236–50 (1976). https://doi.org/10.1137/0713023.
- [Lew82](@cite): J. G. Lewis. *Implementation of the Gibbs–Poole–Stockmeyer and Gibbs–King
    Algorithms*. ACM Transactions on Mathematical Software **8**, 180–89 (1982).
    https://doi.org/10.1145/355993.355998.
"""
struct GibbsPooleStockmeyer <: HeuristicSolver
    node_finder::Function

    function GibbsPooleStockmeyer(node_finder::Function=DEFAULT_NODE_FINDER)
        _assert_valid_node_finder(node_finder)
        return new(node_finder)
    end
end

push!(MatrixBandwidth.ALGORITHMS[:Minimization][:Heuristic], GibbsPooleStockmeyer)

Base.summary(::GibbsPooleStockmeyer) = "Gibbs–Poole–Stockmeyer"

MatrixBandwidth._requires_structural_symmetry(::GibbsPooleStockmeyer) = true

#= We take advantage of the laziness of `Iterators.map` and `Iterators.flatmap` to avoid
allocating `component_orderings` or individual `component[component_ordering]` arrays.
(Indeed, the only allocations performed here are those performed by `connected_components`,
individual `_gps_connected_ordering_reversed` calls, and `collect` at the very end.) =#
function Minimization._minimize_bandwidth_impl(
    A::AbstractMatrix{Bool}, solver::GibbsPooleStockmeyer
)
    node_finder = solver.node_finder
    components = connected_components(A)

    component_orderings = Iterators.map(
        component ->
            _gps_connected_ordering_reversed(view(A, component, component), node_finder),
        components,
    )

    ordering = collect(
        Iterators.flatmap(
            ((component, component_ordering),) -> component[component_ordering],
            zip(components, component_orderings),
        ),
    )
    #= We opt for this unconditional final reversal of the final ordering as opposed to the
    per-component conditional approach detailed in the original Gibbs, Poole, and Stockmeyer
    (1976) paper. This version induces exactly the same bandwidth, missing out on marginal
    profile optimizations in favor of a far simpler implementation. =#
    reverse!(ordering)

    return ordering
end

#= Gibbs–Poole–Stockmeyer searches each connected component independently. Since the usual
Gibbs–Poole–Stockmeyer algorithm reverses the ordering, this "reversed" version does not do
so (i.e., a double reversal, which implies no reversal at all). Instead, the reversal is
carried out in the `_minimize_bandwidth_impl` method instead. =#
function _gps_connected_ordering_reversed(A::AbstractMatrix{Bool}, node_finder::Function)
    u, v = _pseudo_diameter_endpoints(A, node_finder)
    levels_u = _level_structure(A, u)
    levels_v = _level_structure(A, v)
    levels = _combine_level_structures(A, levels_u, levels_v)

    return _number_nodes!(levels, u, v, A)
end

#= Heuristically identify a pair of node that are "farthest" from each other in the
connected component. This is done by generating and iterating over the level structure
generated via a breadth-first search rooted at a node selected by `node_finder`, which
should heuristically identify a node "farthest" from all others in the graph. =#
function _pseudo_diameter_endpoints(A::AbstractMatrix{Bool}, node_finder::Function)
    degrees = vec(sum(A; dims=1))

    u = node_finder(A)
    levels_curr = _level_structure(A, u)
    last_level = levels_curr[end]

    depth_curr = length(levels_curr)
    depths_next = map(node -> length(_level_structure(A, node)), last_level)
    depth_next = maximum(depths_next)

    while depth_next > depth_curr
        max_depth_nodes = last_level[findall(depths_next .== depth_next)]
        u = max_depth_nodes[argmin(degrees[max_depth_nodes])]

        levels_curr = _level_structure(A, u)
        last_level = levels_curr[end]

        depth_curr = depth_next
        depths_next = map(node -> length(_level_structure(A, node)), last_level)
        depth_next = maximum(depths_next)
    end

    widths = Iterators.map(node -> length(_level_structure(A, node)[end]), last_level)
    v = last_level[argmin(widths)]

    return u, v
end

#= Given an adjacency matrix `A`, generate a level partition `L₁, L₂, …, Lₖ` of the graph
represented by `A` rooted at the node `root`. The disjoint union of all levels is precisely
the vertex set of the graph, with each level recursively defined as follows:
    1. `L₁ = {root}`; and
    2. `Lᵢ = {u ∈ N(v) | u ∉ L₁ ∪ L₂ ∪ … ∪ Lᵢ₋₁ and v ∈ Lᵢ₋₁}` for all `i ∈ {2, 3, …, k}`,
        where `N(v)` is the set of neighbors of `v`.
=#
function _level_structure(A::AbstractMatrix{Bool}, root::Int)
    levels = Vector{Int}[]
    level_curr = [root]
    visited = falses(size(A, 1))
    visited[root] = true

    while !isempty(level_curr)
        push!(levels, level_curr)
        level_next = Int[]

        for parent in level_curr
            unvisited = filter(node -> !visited[node], findall(view(A, :, parent)))
            append!(level_next, unvisited)
            visited[unvisited] .= true
        end

        level_curr = level_next
    end

    return levels
end

#= Combine the level structures rooted at `u` and `v` into a single level structure
`L₁, L₂, …, Lₖ` with lower width (i.e., fewer nodes in the largest level) still satisfying
the following properties:
    1. `⋃Lᵢ₌₁ᵏ` is the vertex set of the graph represented by `A`;
    2. for every `u ∈ L₁`, `v ∈ N(u)` implies `v ∈ L₁ ∪ L₂`;
    3. for every `i ∈ {2, 3, …, k - 1}` and `u ∈ Lᵢ`, `v ∈ N(u)` implies
        `v ∈ Lᵢ₋₁ ∪ Lᵢ ∪ Lᵢ₊₁`; and
    4. for every `u ∈ Lₖ`, `v ∈ N(u)` implies `v ∈ Lₖ₋₁ ∪ Lₖ`.
=#
function _combine_level_structures(
    A::AbstractMatrix{Bool}, levels_u::Vector{Vector{Int}}, levels_v::Vector{Vector{Int}}
)
    k = length(levels_u) # Guaranteed to be equal to `length(levels_v)`
    width_u = maximum(Iterators.map(length, levels_u))
    width_v = maximum(Iterators.map(length, levels_v))

    levels = map(_ -> Int[], 1:k)
    level_pairs = Dict(
        node => (findfirst(level -> node in level, levels_v), k - i + 1) for i in 1:k for
        node in levels_u[i]
    )
    assigned = falses(size(A, 1))

    diagonal_pair_nodes = filter(node -> allequal(level_pairs[node]), axes(A, 1))
    A_working = copy(A) # Avoid shared mutability
    A_working[:, diagonal_pair_nodes] .= false
    A_working[diagonal_pair_nodes, :] .= false

    for node in diagonal_pair_nodes
        target_level = level_pairs[node][1]
        push!(levels[target_level], node)
        assigned[node] = true
    end

    components = connected_components(A_working)
    sort!(components; by=length, rev=true)

    for component in connected_components(A_working)
        sizes_curr = length.(levels)
        pot_sizes_u = copy(sizes_curr)
        pot_sizes_v = copy(sizes_curr)

        component_level_pairs = map(node -> level_pairs[node], component)
        component_level_pairs_u = first.(component_level_pairs)
        component_level_pairs_v = last.(component_level_pairs)
        pot_sizes_u[component_level_pairs_u] .+= 1
        pot_sizes_v[component_level_pairs_v] .+= 1

        if !isempty(component_level_pairs_u)
            max_size_u = maximum(view(pot_sizes_u, component_level_pairs_u))
        else
            max_size_u = 0
        end

        if !isempty(component_level_pairs_v)
            max_size_v = maximum(view(pot_sizes_v, component_level_pairs_v))
        else
            max_size_v = 0
        end

        if max_size_u < max_size_v || (max_size_u == max_size_v && width_u <= width_v)
            target_idx = 1
        else
            target_idx = 2
        end

        for node in filter(node -> !assigned[node], component)
            target_level = level_pairs[node][target_idx]
            push!(levels[target_level], node)
            assigned[node] = true
        end
    end

    return levels
end

#= Number the nodes in the level structure `levels` in such a way that is likely to minimize
the bandwidth of `A`. =#
function _number_nodes!(
    levels::Vector{Vector{Int}}, u::Int, v::Int, A::AbstractMatrix{Bool}
)
    n = size(A, 1)

    ordering = Vector{Int}(undef, size(A, 1))
    degrees = vec(sum(A; dims=1))

    if degrees[u] < degrees[v]
        reverse!(levels)
        start = u
    else
        start = v
    end

    ordering[1] = start
    num_placed = 1

    for (i, level) in enumerate(levels)
        placed_nodes = view(ordering, 1:num_placed)
        unvisited = Set(level)

        if i == 1
            delete!(unvisited, start) # `start` is guaranteed to be in `levels[1]`
            #= Just a placeholder value for the sorting step, since `nums_dangling` is not
            actually needed for the first level. =#
            nums_dangling = falses(n)
        else
            nums_dangling = sum(view(A, placed_nodes, :); dims=1)
        end

        while !isempty(unvisited)
            #= While placing the nodes of the first level, we update `placed_nodes`
            simultaneously as we go. =#
            if i == 1
                placed_nodes = view(ordering, 1:num_placed)
            end

            lowest_prev_with_unvisited_neighbor = 0
            min_ordering_num = n + 1
            res = iterate(placed_nodes)
            idx = 1

            while (!isnothing(res) && idx < min_ordering_num)
                (node, state) = res

                if !isempty(intersect(findall(view(A, :, node)), unvisited))
                    lowest_prev_with_unvisited_neighbor = node
                    min_ordering_num = idx
                end

                res = iterate(placed_nodes, state)
                idx += 1
            end

            if lowest_prev_with_unvisited_neighbor == 0
                min_degree_unvisited = reduce(
                    (x, y) -> (degrees[x], x) < (degrees[y], y) ? x : y, unvisited
                )
                ordering[num_placed += 1] = min_degree_unvisited
                delete!(unvisited, min_degree_unvisited)
            else
                unvisited_adj = intersect(
                    findall(view(A, :, lowest_prev_with_unvisited_neighbor)), unvisited
                )
                sort!(unvisited_adj; by=node -> (-nums_dangling[node], degrees[node]))

                copyto!(ordering, num_placed + 1, unvisited_adj)
                num_placed += length(unvisited_adj)
                filter!(!in(unvisited_adj), unvisited)
            end
        end
    end

    return ordering
end
