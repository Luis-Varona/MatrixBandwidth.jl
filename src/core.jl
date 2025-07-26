# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    bandwidth(A) -> Int

Compute the bandwidth of `A` before any permutation of its rows and columns.

The *bandwidth* of an ``n×n`` matrix ``A`` is the minimum non-negative integer
``k ∈ [0, n - 1]`` such that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A``
has bandwidth *at most* ``k`` if all entries above the ``k``ᵗʰ superdiagonal and below the
``k``ᵗʰ subdiagonal are zero, and ``A`` has bandwidth *at least* ``k`` if there exists any
nonzero entry in the ``k``ᵗʰ superdiagonal or subdiagonal.

In contrast to [`minimize_bandwidth`](@ref), this function does not attempt to minimize the
bandwidth of `A` by permuting its rows and columns—it simply computes its bandwidth as is.

# Arguments
- `A::AbstractMatrix{<:Number}`: the (square) matrix whose bandwidth is computed.

# Returns
- `::Int`: the bandwidth of `A`.

# Performance
Given an ``n×n`` input matrix ``A``, this relatively simple algorithm runs in ``O(n²)``
time.

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

    indices = findall(!iszero, A)

    if isempty(indices)
        band = 0 # Avoid reducing over an empty collection with `maximum``
    else
        band = maximum(Iterators.map(idx -> abs(idx[1] - idx[2]), indices))
    end

    return band
end

"""
    profile(A) -> Int

Compute the profile of `A` before any permutation of its rows and columns.

The *profile* of a structurally symmetric ``n×n`` matrix ``A`` is traditionally defined as
the sum of the distances from each diagonal entry to the leftmost nonzero entry in that
row—in other words, ``∑ᵢ₌₁ⁿ (i - fᵢ)``, where each ``fᵢ`` is the smallest index such that
``A[i, fᵢ] ≠ 0`` [Maf14; pp. 187-88](@cite). Generalizing this property to all square
matrices, we define the *column profile* of a matrix to be the sum of the distances from
each diagonal entry to the farthest (not necessarily topmost) nonzero entry in that column
and the *row profile* to be the sum of the distances from each diagonal entry to the
farthest (not necessarily leftmost) nonzero entry in that row. (Note that both of these
properties are equal to traditional matrix profile for structurally symmetric matrices.)

One of the most common contexts in which matrix profile is relevant is sparse matrix
storage, where lower-profile matrices occupy less space in memory [Maf14; p.188](@cite).
Since Julia's `SparseArrays` package defaults to compressed sparse column storage over
compressed sparse row, we therefore compute column profile by default unless the dimension
is otherwise specified.

# Arguments
- `A::AbstractMatrix{<:Number}`: the (square) matrix whose profile is computed.

# Keyword Arguments
- `dim::Symbol=:col`: the dimension along which the profile is computed; must be either
    `:col` (the default) or `:row`.

# Returns
- `::Int`: the profile of `A` along the specified dimension.

# Performance
Given an ``n×n`` input matrix ``A``, this relatively simple algorithm runs in ``O(n²)``
time.

# Examples
`profile` computes the column profile of a matrix by default:
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(2287);

julia> (n, p) = (25, 0.05);

julia> A = sprand(n, n, p)
25×25 SparseMatrixCSC{Float64, Int64} with 29 stored entries:
⎡⠀⠀⠀⠀⠀⠀⠐⠀⠒⠀⡀⠀⠀⎤
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠐⠂⎥
⎢⠀⠀⠀⢀⠌⠀⠀⠀⢀⠈⠀⠀⠀⎥
⎢⠠⢄⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠂⡀⠀⠀⠌⠀⠈⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠢⡀⢄⡈⠀⠀⠀⎥
⎣⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⎦

julia> profile(A)
211
```

The dimension (`:row` or `:col`) can also be explicitly specified:
```jldoctest
julia> using Random, SparseArrays

julia> Random.seed!(3647);

julia> (n, p) = (25, 0.05);

julia> A = sprand(n, n, p)
25×25 SparseMatrixCSC{Float64, Int64} with 31 stored entries:
⎡⠄⠀⠀⠀⠀⠀⠀⠘⠀⠀⠀⠁⠀⎤
⎢⠄⢀⠀⠀⠁⠀⠀⠀⠀⢀⠀⠀⠀⎥
⎢⠀⢀⡂⠀⠀⠀⠀⠀⠀⠀⠀⡄⠀⎥
⎢⠀⠀⠀⠀⠀⡀⠂⠀⠀⠀⠀⠀⠀⎥
⎢⠁⠀⠁⠀⠁⠀⠀⠁⠄⢀⠈⠀⠀⎥
⎢⠂⠐⠐⠐⠠⠀⠄⠀⠀⠀⠠⣀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎦

julia> profile(A, dim=:row)
147

julia> profile(A, dim=:col)
175
```
"""
function profile(A::AbstractMatrix{<:Number}; dim::Symbol=:col)
    if dim == :col
        slice_iter = eachcol(A)
    elseif dim == :row
        slice_iter = eachrow(A)
    else
        throw(ArgumentError("Expected `dim` to be `:col` or `:row`, got `$dim`"))
    end

    m, n = size(A)

    if m != n
        throw(RectangularMatrixError(A))
    end

    prof = 0

    for (i, slice) in enumerate(slice_iter)
        j = max(i - 1, n - i + 1)
        found_nonzero = false

        while (!found_nonzero && j > 0)
            if i - j > 0 && slice[i - j] != 0
                prof += j
                found_nonzero = true
            elseif i + j <= n && slice[i + j] != 0
                prof += j
                found_nonzero = true
            end

            j -= 1
        end
    end

    return prof
end

"""
    bandwidth_lower_bound(A) -> Int

Compute a lower bound on the bandwidth of `A` using [CSG05; pp. 359--60](@cite)'s results.

`A` is assumed to be structurally symmetric, since the bound from
[CSG05; pp.359--60](@cite) was discovered in the context of undirected graphs (whose
adjacency matrices are symmetric).

The *bandwidth* of an ``n×n`` matrix ``A`` is the minimum non-negative integer
``k ∈ [0, n - 1]`` such that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A``
has bandwidth *at most* ``k`` if all entries above the ``k``ᵗʰ superdiagonal and below the
``k``ᵗʰ subdiagonal are zero, and ``A`` has bandwidth *at least* ``k`` if there exists any
nonzero entry in the ``k``ᵗʰ superdiagonal or subdiagonal.

In contrast to [`minimize_bandwidth`](@ref), this function does not attempt to truly
minimize the bandwidth of `A`—it simply returns a lower bound on its bandwidth up to
symmetric permutation of its rows and columns. This bound is not generally tight, but it
indeed matches the true minimum in many non-trivial cases and is easily computable in
``O(n³)`` time (dominated by the Floyd–Warshall algorithm call; the core logic itself runs
in ``O(n²)`` time).

# Arguments
- `A::AbstractMatrix{<:Number}`: the (square) matrix on whose bandwidth a lower bound is to
    be computed. `A` must be structurally symmetric (i.e., `A[i, j]` must be nonzero if
    and only if `A[j, i]` is nonzero for ``1 ≤ i, j ≤ n``).

# Returns
- `::Int`: a lower bound on the bandwidth of `A`. (This bound is tight in many non-trivial
    cases but not universally so.)

# Examples
The function correctly computes a bound less than (or equal to) the true minimum bandwidth
of a matrix up to symmetric permutation:
```jldoctest
julia> using Random, SparseArrays, Combinatorics

julia> Random.seed!(21);

julia> (n, p) = (9, 0.4);

julia> A = sprand(n, n, p);

julia> A = A + A' # Ensure structural symmetry;

julia> minimize_bandwidth(A, Minimization.BruteForceSearch())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Brute-force search
 * Approach: exact
 * Minimum Bandwidth: 5
 * Original Bandwidth: 8
 * Matrix Size: 9×9

julia> bandwidth_lower_bound(A) # Always less than or equal to the true minimum bandwidth
4
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

    #= The bandwidth is trivially zero, so we return early. This also prevents `gamma` from
    erroneously never being updated to 0.0, since `finite_nonzero_dists` would always be
    empty if we were to skip this check. =#
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
        # Exclude not only unreachable nodes but also self-distances (always 0)
        finite_nonzero_dists = Int.(filter(k -> isfinite(k) && k != 0.0, dists))

        # Only non-isolated nodes are used to update `alpha` and `gamma`
        if !isempty(finite_nonzero_dists)
            max_dist = maximum(finite_nonzero_dists)
            alpha_cand = 0
            gamma_cand = 0

            #= Compute the `k`-hop neighborhood sizes for each distance `k ≥ 1`. Every node
            is a 0-hop neighbor of itself, so we add 1 after taking the cumulative sum
            (since we filtered out self-distances in `finite_nonzero_dists`). =#
            k_hop_neighborhood_sizes = zeros(Int, max_dist)
            foreach(k -> k_hop_neighborhood_sizes[k] += 1, finite_nonzero_dists)
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
