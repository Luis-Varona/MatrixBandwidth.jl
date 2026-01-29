# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    random_banded_matrix(n, k; p=0.5, rng=default_rng()) -> Matrix{Float64}

Generate a random `n×n` structurally symmetric `k`-banded matrix with band density `≈ p`.

By definition of structural symmetry, the ``(i, j)``-th entry of the matrix is nonzero if
and only if the ``(j, i)``-th entry is nonzero as well. All entries from this matrix are
from the interval `[0, 1]`. Entries up to the `k`ᵗʰ superdiagonal and down to the `k`ᵗʰ
subdiagonal are nonzero with probability `p`.

It is also guaranteed that each of these bands (besides the main diagonal) has at least one
nonzero entry (even when `p` is very small), thus ensuring that the matrix has bandwidth
precisely `k` before any reordering. (There may, however, still exist a symmetric
permutation inducing a minimum bandwidth less than `k`, especially for small values of `p`.)

# Arguments
- `n::Integer`: the order of the matrix to generate. Must be positive.
- `k::Integer`: the desired matrix bandwidth. Must satisfy `0 ≤ k < n`.

# Keyword Arguments
- `p::Real=0.5`: the band density. Must satisfy `0 < p ≤ 1`. Defaults to `0.5`.
- `rng::AbstractRNG=Random.default_rng()`: the random number generator to use. Defaults to
    `Random.default_rng()`.

# Returns
- `::Matrix{Float64}`: a random `n×n` matrix with bandwidth exactly `k` and sparse bands
    with density `p`.

# Examples
Generate a ``6×6`` matrix with bandwidth ``1`` and the maximum number of nonzero entries:
```jldoctest
julia> using Random

julia> A = MatrixBandwidth.random_banded_matrix(6, 1; p=1, rng=MersenneTwister(1228))
6×6 Matrix{Float64}:
 0.310239  0.346413  0.0       0.0        0.0       0.0
 0.509981  0.917073  0.390771  0.0        0.0       0.0
 0.0       0.760045  0.808396  0.0195686  0.0       0.0
 0.0       0.0       0.222338  0.853164   0.806888  0.0
 0.0       0.0       0.0       0.421603   0.132165  0.805813
 0.0       0.0       0.0       0.0        0.305339  0.0799183

julia> bandwidth(A)
1
```

Generate a ``7×7`` matrix with bandwidth ``3`` and band density ``0.3``:
```jldoctest
julia> using Random

julia> A = MatrixBandwidth.random_banded_matrix(7, 3; p=0.3, rng=MersenneTwister(0402))
7×7 Matrix{Float64}:
 0.0       0.132699  0.0       0.0       0.0  0.0       0.0
 0.869352  0.0       0.324319  0.926496  0.0  0.0       0.0
 0.0       0.891878  0.0       0.658102  0.0  0.0       0.0
 0.0       0.88859   0.399559  0.0       0.0  0.284285  0.703377
 0.0       0.0       0.0       0.0       0.0  0.0       0.0
 0.0       0.0       0.0       0.489594  0.0  0.0       0.393573
 0.0       0.0       0.0       0.412412  0.0  0.47063   0.0

julia> bandwidth(A)
3
```

Generate an ``8×8`` diagonal (bandwidth ``0``) matrix with default band density (``0.5``):
```jldoctest
julia> using Random

julia> A = MatrixBandwidth.random_banded_matrix(8, 0; rng=MersenneTwister(0102))
8×8 Matrix{Float64}:
 0.0  0.0        0.0       0.0       0.0  0.0      0.0  0.0
 0.0  0.0762399  0.0       0.0       0.0  0.0      0.0  0.0
 0.0  0.0        0.373113  0.0       0.0  0.0      0.0  0.0
 0.0  0.0        0.0       0.726309  0.0  0.0      0.0  0.0
 0.0  0.0        0.0       0.0       0.0  0.0      0.0  0.0
 0.0  0.0        0.0       0.0       0.0  0.41974  0.0  0.0
 0.0  0.0        0.0       0.0       0.0  0.0      0.0  0.0
 0.0  0.0        0.0       0.0       0.0  0.0      0.0  0.293132

julia> bandwidth(A)
0
```

# Notes
Users of the [`MatrixBandwidth`](@ref) package may find this function useful when generating
random test data for whatever frameworks, algorithms, etc. they are implementing.
"""
function random_banded_matrix(
    n::Integer, k::Integer; p::Real=0.5, rng::AbstractRNG=Random.default_rng()
)
    if n <= 0
        throw(ArgumentError("Matrix order must be positive, got $n"))
    end

    if k < 0
        throw(ArgumentError("Matrix bandwidth must be non-negative, got $k"))
    end

    if n <= k
        throw(ArgumentError("Matrix order must be greater than bandwidth, got $n and $k"))
    end

    if !(0 < p <= 1)
        throw(ArgumentError("Band density must be in (0, 1], got $p"))
    end

    A = zeros(Float64, n, n)

    #= The main diagonal is handled separately to halve the number of operations and avoid
    the need to ensure at least one nonzero entry (as this would not affect bandwidth). =#
    for i in 1:n
        if rand(rng) < p
            A[i, i] = rand(rng)
        end
    end

    #= Conditionally add entries to the superdiagonals and subdiagonals (up to the `k`ᵗʰ)
    based on the band density `p`. =#
    for d in 1:k
        feasible_indices = 1:(n - d)

        for i in feasible_indices
            if rand(rng) < p
                A[i, i + d] = rand(rng)
                A[i + d, i] = rand(rng)
            end
        end

        #= Ensure that the `d`ᵗʰ superdiagonal and subdiagonal each have at least one
        nonzero entry for `1 ≤ d ≤ k`, making a `bandwidth < k` symmetric permutation less
        likely. (Given structural symmetry, we need only check for empty superdiagonals.) =#
        if all(iszero, Iterators.map(i -> A[i, i + d], feasible_indices))
            i = rand(rng, feasible_indices)
            A[i, i + d] = rand(rng)
            A[i + d, i] = rand(rng)
        end
    end

    return A
end

"""
    connected_components(A) -> Vector{Vector{Int}}

Find the indices of all connected components of the graph whose adjacency matrix is `A`.

`A` is assumed to be symmetric, representing an undirected graph.

# Arguments
- `A::AbstractMatrix{Bool}`: the adjacency matrix of the graph. Must be symmetric.

# Returns
- `::Vector{Vector{Int}}`: a vector of vectors, where each element is a vector of node
    indices belonging to a connected component.

# Examples
```jldoctest
julia> using Graphs

julia> g = complement(complete_multipartite_graph([3, 4, 2]))
{9, 10} undirected simple Int64 graph

julia> A = Bool.(adjacency_matrix(g))
9×9 SparseArrays.SparseMatrixCSC{Bool, Int64} with 20 stored entries:
 ⋅  1  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 1  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 1  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  1  1  1  ⋅  ⋅
 ⋅  ⋅  ⋅  1  ⋅  1  1  ⋅  ⋅
 ⋅  ⋅  ⋅  1  1  ⋅  1  ⋅  ⋅
 ⋅  ⋅  ⋅  1  1  1  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅

julia> MatrixBandwidth.connected_components(A)
3-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [4, 5, 6, 7]
 [8, 9]
```
"""
function connected_components(A::AbstractMatrix{Bool})
    n = size(A, 1)
    visited = falses(n)
    queue = Queue{Int}()
    components = Vector{Int}[]

    for i in 1:n
        if !visited[i]
            visited[i] = true
            push!(queue, i)
            component = Int[]

            while !isempty(queue)
                u = popfirst!(queue)
                push!(component, u)

                for v in findall(view(A, :, u))
                    if !visited[v]
                        visited[v] = true
                        push!(queue, v)
                    end
                end
            end

            push!(components, component)
        end
    end

    return components
end

"""
    floyd_warshall_shortest_paths(A) -> Matrix{Float64}

Compute a distance matrix from the adjacency matrix `A` using the Floyd–Warshall algorithm.

Relatively isolated pairs of nodes (those unreachable from each other) are assigned
distances of `Inf`.

`A` is assumed to be symmetric with an all-false diagonal, representing a simple graph.

# Arguments
- `A::AbstractMatrix{Bool}`: the adjacency matrix of the graph. Must be symmetric with an
    all-false diagonal.

# Returns
- `::Matrix{Float64}`: an `n×n` matrix `D`, where `D[i, j]` is the length of the shortest
    path between nodes `i` and `j`, or `Inf` if no such path exists.

# Performance
Given an ``n×n`` input matrix, the Floyd–Warshall algorithm runs in ``O(n³)`` time, given
that each level in the triple-nested loop iterates over ``O(n)`` entries.

# Examples
Floyd–Warshall finds the shortest distances between all pairs of nodes in a connected graph:
```jldoctest
julia> using Graphs

julia> g = ladder_graph(5)
{10, 13} undirected simple Int64 graph

julia> A = Bool.(adjacency_matrix(g))
10×10 SparseArrays.SparseMatrixCSC{Bool, Int64} with 26 stored entries:
 ⋅  1  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅
 1  ⋅  1  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅
 ⋅  1  ⋅  1  ⋅  ⋅  ⋅  1  ⋅  ⋅
 ⋅  ⋅  1  ⋅  1  ⋅  ⋅  ⋅  1  ⋅
 ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  1
 1  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅
 ⋅  1  ⋅  ⋅  ⋅  1  ⋅  1  ⋅  ⋅
 ⋅  ⋅  1  ⋅  ⋅  ⋅  1  ⋅  1  ⋅
 ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  1  ⋅  1
 ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  1  ⋅

julia> Int.(MatrixBandwidth.floyd_warshall_shortest_paths(A))
10×10 Matrix{Int64}:
 0  1  2  3  4  1  2  3  4  5
 1  0  1  2  3  2  1  2  3  4
 2  1  0  1  2  3  2  1  2  3
 3  2  1  0  1  4  3  2  1  2
 4  3  2  1  0  5  4  3  2  1
 1  2  3  4  5  0  1  2  3  4
 2  1  2  3  4  1  0  1  2  3
 3  2  1  2  3  2  1  0  1  2
 4  3  2  1  2  3  2  1  0  1
 5  4  3  2  1  4  3  2  1  0
```

Floyd–Warshall assigns `Inf` to pairs of nodes in different connected components:
```jldoctest
julia> using Graphs

julia> g = complement(wheel_graph(8))
{8, 14} undirected simple Int64 graph

julia> A = Bool.(adjacency_matrix(g))
8×8 SparseArrays.SparseMatrixCSC{Bool, Int64} with 28 stored entries:
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  1  1  1  1  ⋅
 ⋅  ⋅  ⋅  ⋅  1  1  1  1
 ⋅  1  ⋅  ⋅  ⋅  1  1  1
 ⋅  1  1  ⋅  ⋅  ⋅  1  1
 ⋅  1  1  1  ⋅  ⋅  ⋅  1
 ⋅  1  1  1  1  ⋅  ⋅  ⋅
 ⋅  ⋅  1  1  1  1  ⋅  ⋅

julia> MatrixBandwidth.floyd_warshall_shortest_paths(A)
8×8 Matrix{Float64}:
  0.0  Inf   Inf   Inf   Inf   Inf   Inf   Inf
 Inf    0.0   2.0   1.0   1.0   1.0   1.0   2.0
 Inf    2.0   0.0   2.0   1.0   1.0   1.0   1.0
 Inf    1.0   2.0   0.0   2.0   1.0   1.0   1.0
 Inf    1.0   1.0   2.0   0.0   2.0   1.0   1.0
 Inf    1.0   1.0   1.0   2.0   0.0   2.0   1.0
 Inf    1.0   1.0   1.0   1.0   2.0   0.0   2.0
 Inf    2.0   1.0   1.0   1.0   1.0   2.0   0.0
```
"""
function floyd_warshall_shortest_paths(A::AbstractMatrix{Bool})
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

"""
    is_structurally_symmetric(A) -> Bool

Check whether `A[i, j]` is nonzero if and only if `A[j, i]` is nonzero for all `i` and `j`.

# Arguments
- `A::AbstractMatrix{<:Number}`: the matrix to check for structural symmetry.

# Returns
- `::Bool`: whether `A` is structurally symmetric.

# Examples
```jldoctest
julia> A = [4 0 9 -2; 0 0 1 0; 3 -1 5 0; 4 0 0 3]
4×4 Matrix{Int64}:
 4   0  9  -2
 0   0  1   0
 3  -1  5   0
 4   0  0   3

julia> MatrixBandwidth.is_structurally_symmetric(A)
true

julia> B = [1.12 2.36 0.00; 5.99 0.0 0.0; 0.0 3.1 -7.49]
3×3 Matrix{Float64}:
 1.12  2.36   0.0
 5.99  0.0    0.0
 0.0   3.1   -7.49

julia> MatrixBandwidth.is_structurally_symmetric(B)
false
```

# Notes
Instead of transposing `A` and allocating a new matrix, it suffices to iterate over all
opposing pairs of off-diagonal entries.
"""
function is_structurally_symmetric(A::AbstractMatrix{<:Number})
    (m, n) = size(A)
    return m == n &&
           all((A[i, j] != 0) == (A[j, i] != 0) for i in 1:(n - 1) for j in (i + 1):n)
end

"""
    offdiag_nz_support(A) -> AbstractMatrix{Bool}

Convert `A` to a boolean matrix and set all its diagonal entries to `false`.

# Arguments
- `A::AbstractMatrix{<:Number}`: the matrix to convert.

# Returns
- `::AbstractMatrix{Bool}`: the nonzero support of `A`, with all diagonal entries set to
    `false`.

# Examples
```jldoctest
julia> A = [0 2 0 7; 0 -8 0 3; -1 9 0 0; 0 0 0 5]
4×4 Matrix{Int64}:
  0   2  0  7
  0  -8  0  3
 -1   9  0  0
  0   0  0  5

julia> MatrixBandwidth.offdiag_nz_support(A)
4×4 BitMatrix:
 0  1  0  1
 0  0  0  1
 1  1  0  0
 0  0  0  0
```

# Notes
In the context of matrix bandwidth reduction algorithms (which are only concerned with the
nonzero support of the input matrix), this improves performance via cache optimizations,
availability of bitwise operations, etc.
"""
function offdiag_nz_support(A::AbstractMatrix{<:Number})
    A_bool = (!iszero).(A)
    foreach(i -> A_bool[i, i] = false, axes(A_bool, 1))
    return A_bool
end

"""
    find_direct_subtype(abstracttype, subtype) -> Type

Identify the highest supertype of `subtype` that is also a subtype of `abstracttype`.

# Arguments
- `abstracttype::Type`: an abstract type.
- `subtype::Type`: a subtype of `abstracttype`.

# Returns
- `::Type`: the direct subtype of `abstracttype` that is a supertype of `subtype`.

# Examples
```jldoctest
julia> abstract type Parent end

julia> abstract type Child1 <: Parent end

julia> abstract type Grandchild1 <: Child1 end

julia> struct Grandchild2 <: Child1 end

julia> abstract type Child2 <: Parent end

julia> struct Child3 <: Parent end

julia> MatrixBandwidth.find_direct_subtype(Parent, Child1)
Child1

julia> MatrixBandwidth.find_direct_subtype(Parent, Grandchild1)
Child1

julia> MatrixBandwidth.find_direct_subtype(Parent, Grandchild2)
Child1

julia> MatrixBandwidth.find_direct_subtype(Parent, Child2)
Child2

julia> MatrixBandwidth.find_direct_subtype(Parent, Child3)
Child3
```
"""
function find_direct_subtype(abstracttype::Type, subtype::Type)
    if !isabstracttype(abstracttype)
        throw(ArgumentError("Expected an abstract type, got $abstracttype"))
    end

    if !(subtype <: abstracttype)
        throw(ArgumentError("Expected a subtype of $abstracttype, got $subtype"))
    end

    parent = supertype(subtype)

    if parent === abstracttype
        child = subtype
    else
        child = find_direct_subtype(abstracttype, parent)
    end

    return child
end
