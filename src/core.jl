# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    bandwidth(A) -> Int

Compute the bandwidth of `A` before any permutation of its rows and columns.

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

In contrast to [`minimize_bandwidth`](@ref), this function does not attempt to minimize the
bandwidth of `A` by permuting its rows and columns—it simply computes its bandwidth as is.

# Arguments
- `A::AbstractMatrix{<:Number}`: the (square) matrix whose bandwidth is computed.

# Returns
- `::Int`: the bandwidth of `A`.

# Performance
Given an ``n×n`` input matrix ``A``, this relatively simple algorithm runs in ``O(n²)`` time.

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
 1  0  0  0  0  0  0  0
 0  1  0  1  0  0  0  0
 0  0  0  1  1  0  0  0
 0  1  1  1  0  1  0  0
 0  0  1  0  0  0  0  0
 0  0  0  1  0  0  0  0
 0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0

julia> bandwidth(A)
2

julia> A_shuffled = A[perm, perm]
8×8 BitMatrix:
 0  0  0  0  0  0  1  0
 0  0  0  0  0  0  0  0
 0  0  0  1  0  0  0  0
 0  0  1  0  0  0  1  0
 0  0  0  0  1  0  0  0
 0  0  0  0  0  1  1  0
 1  0  0  1  0  1  1  0
 0  0  0  0  0  0  0  0

julia> bandwidth(A_shuffled)
6
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
function bandwidth(A::AbstractMatrix{<:Number})
    if !allequal(size(A))
        throw(RectangularMatrixError(A))
    end

    #= We are only concerned with which (off-diagonal) entries are nonzero, not the actual
    values. We also set every diagonal entry to `false` for consistency with any algorithms
    that assume an adjacency matrix structure. =#
    A_bool = _offdiag_nonzero_support(A)

    indices = findall(A_bool)

    if isempty(indices)
        band = 0
    else
        band = maximum(abs(index[1] - index[2]) for index in indices)
    end

    return band
end

"""
    bandwidth_lower_bound(A) -> Int

Compute a lower bound on the bandwidth of `A` using [CSG05](@cite)'s results.

The nonzero support of `A` is assumed to be symmetric, since [CSG05](@cite)'s bound was
discovered in the context of undirected graphs (whose adjacency matrices are symmetric).

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

In contrast to [`minimize_bandwidth`](@ref), this function does not attempt to truly
minimize the bandwidth of `A`—it simply returns a lower bound on its bandwidth up to
symmetric permutation of its rows and columns. This bound is not tight, but it is easily
computable in `O(n³)` time, dominated by the Floyd–Warshall algorithm call. (The core logic
here runs in `O(n²)` time.)

# Arguments
- `A::AbstractMatrix{<:Number}`: the (square) matrix on whose bandwidth a lower bound is to
    be computed. `A` must have a symmetric nonzero support (i.e., `A[i, j]` is nonzero if
    and only if `A[j, i]` is nonzero).

# Returns
- `::Int`: a lower bound on the bandwidth of `A`. (This bound is not tight.)

# Examples
TODO: Write here

# Notes
Some texts define matrix bandwidth to be the minimum non-negative integer ``k`` such that
``A[i, j] = 0`` whenever ``|i - j| ≥ k`` instead, particularly in more mathematically-minded
communities. Effectively, this definition treats diagonal matrices as bandwidth ``1``,
tridiagonal matrices as bandwidth ``2``, and so on. Our definition, on the other hand, is
more common in computer science contexts, treating diagonal matrices as bandwidth ``0`` and
tridiagonal matrices as bandwidth ``1``. (Both definitions, however, agree that the
bandwidth of an empty matrix is simply ``0``.)
"""
function bandwidth_lower_bound(A::AbstractMatrix{<:Number})
    if !allequal(size(A))
        throw(RectangularMatrixError(A))
    end

    if !_is_structurally_symmetric(A)
        throw(
            DomainError(
                A,
                "Bandwidth lower bound is only computable for structurally symmetric matrices",
            ),
        )
    end

    #= We are only concerned with which (off-diagonal) entries are nonzero, not the actual
    values. We also set every diagonal entry to `false` for consistency with any algorithms
    that assume an adjacency matrix structure. =#
    A_bool = _offdiag_nonzero_support(A)

    #= The bandwidth is trivially zero, so we return early. This also anticipates `gamma`
    never being updated to 0.0, since in that case, `length(finite_dists)` will always be
    precisely 1 (due to the diagonal entries of `dist_matrix`). =#
    if iszero(A_bool)
        return 0
    end

    n = size(A_bool, 1)

    # The bandwidth is trivially `n - 1`, so we return early
    if all(A_bool[i, j] for i in 1:(n - 1) for j in (i + 1):n)
        return n - 1
    end

    # If `A_bool` is not symmetric, an error is thrown here
    dist_matrix = _floyd_warshall_shortest_paths(A_bool)
    alpha = 0 # The minimum possible bandwidth is 0
    gamma = n - 1 # The maximum possible bandwidth is `n - 1`

    for dists in eachcol(dist_matrix)
        finite_dists = Int.(filter(isfinite, dists)) # Exclude unreachable nodes

        #= Since the distance from a node to itself is always 0, `finite_dists` will always
        contain at least one element. Only non-isolated nodes are used to update `alpha` and
        `gamma`, so we check that `finite_dists` contains more than just this element. =#
        if length(finite_dists) > 1
            max_dist = maximum(finite_dists)
            alpha_cand = 0
            gamma_cand = 0

            #= Compute the `k`-hop neighborhood sizes for each distance `k ≥ 1`. Every node
            is a 0-hop neighbor of itself, so we add 1 after taking the cumulative sum. =#
            k_hop_neighborhood_sizes = zeros(Int, max_dist)
            foreach(k -> k_hop_neighborhood_sizes[k] += 1, filter(!iszero, finite_dists))
            k_hop_neighborhood_sizes .= cumsum(k_hop_neighborhood_sizes) .+ 1

            for (k, num_k_hop_neighbors) in enumerate(k_hop_neighborhood_sizes)
                #= Every node is a 0-hop neighbor of itself, so subtract 1. (We could have
                just avoided adding 1 when constructing `k_hop_neighborhood_sizes`, but this
                implementation aligns better with the graph-theoretic definition. =#
                alpha_cand = max(alpha_cand, cld(num_k_hop_neighbors - 1, 2k))
                gamma_cand = max(gamma_cand, cld(num_k_hop_neighbors - 1, k))
            end

            alpha = max(alpha, alpha_cand)
            gamma = min(gamma, gamma_cand)
        end
    end

    return max(alpha, gamma) # Take the better lower bound on the matrix bandwidth
end

#= Compute a distance matrix from an `n×n` adjacency matrix `A` in `O(n³)` time. Relatively
isolated pairs of nodes (unreachable from each other) are assigned a distance of `Inf`. `A`
is assumed to be symmatric with an all-false diagonal. =#
function _floyd_warshall_shortest_paths(A::AbstractMatrix{Bool})
    n = size(A, 1)

    D = Matrix{Float64}(undef, n, n)
    foreach(i -> D[i, i] = 0.0, 1:n)

    for i in 1:(n - 1), j in (i + 1):n
        if A[i, j]
            D[i, j] = D[j, i] = 1.0
        else
            D[i, j] = D[j, i] = Inf
        end
    end

    for k in 1:n, i in 1:(n - 1), j in (i + 1):n
        if D[i, j] > D[i, k] + D[k, j]
            D[i, j] = D[j, i] = D[i, k] + D[k, j]
        end
    end

    return D
end
