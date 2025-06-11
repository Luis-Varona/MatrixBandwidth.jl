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
produced by application of Cuthill–McKee; [Geo71; pp. 114--15](@cite) found that this tends
to produce an even better ordering.

We also extend the algorithm to work more generally when ``A`` is not symmetric by applying
it to ``A + Aᵀ`` instead, as suggested in [RS06; p. 808](@cite). This approach still tends
to produce a fairly good ordering, but it is not guaranteed to be as optimal as directly
applying reverse Cuthill–McKee to a symmetric input.

# Fields
- `node_selector::Function`: a function that selects a node from some connected component of
    the input matrix from which to start the breadth-first search. If no custom heuristic is
    specified, this field defaults to [`pseudo_peripheral_node`](@ref), which picks a node
    "farthest" from the others in the component (not necessarily the lowest-degree node).

# Examples
Reverse Cuthill–McKee finds an optimal ordering for an asymmetric ``45×45`` matrix with
bandwidth ``4`` whose rows and columns have been shuffled:
```jldoctest
julia> using Random

julia> Random.seed!(87);

julia> (n, k) = (45, 4);

julia> perm = randperm(n);

julia> A = MatrixBandwidth.random_sparse_banded_matrix(n, k);

julia> A_shuffled = A[perm, perm];

julia> res = minimize_bandwidth(A_shuffled, ReverseCuthillMcKee());

julia> iszero.(A_shuffled) == iszero.(A_shuffled') # Works even for asymmetric matrices
false

julia> bandwidth_unpermuted(A)
4

julia> bandwidth_unpermuted(A_shuffled)
43

julia> res.bandwidth # The true minimum bandwidth
4
```

Reverse Cuthill–McKee finds a near-optimal ordering for an asymmetric ``250×250`` matrix
with bandwidth ``14`` whose rows and columns have been shuffled:
```jldoctest
julia> using Random

julia> Random.seed!(5747);

julia> (n, k) = (250, 14);

julia> perm = randperm(n);

julia> A = MatrixBandwidth.random_sparse_banded_matrix(n, k);

julia> A_shuffled = A[perm, perm];

julia> res = minimize_bandwidth(A_shuffled, ReverseCuthillMcKee());

julia> iszero.(A_shuffled) == iszero.(A_shuffled') # Works even for asymmetric matrices
false

julia> bandwidth_unpermuted(A)
14

julia> bandwidth_unpermuted(A_shuffled)
245

julia> res.bandwidth # Close to the true minimum
17
```

# Notes
See [`CuthillMcKee`](@ref) and the associated method of `_bool_minimal_band_ordering` for
our implementation of the original Cuthill–McKee algorithm. (Indeed, the reverse
Cuthill–McKee method of `_bool_minimal_band_ordering` is merely a wrapper around the
Cuthill–McKee method.)

Note also that the `node_selector` field must be of the form
`(A::AbstractMatrix{Bool}) -> Integer` (i.e., it must take in an boolean matrix and return
an integer). If this is not the case, an `ArgumentError` is thrown upon construction.

See also the documentation for supertypes [`HeuristicSolver`](@ref) and
[`AbstractSolver`](@ref).
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
