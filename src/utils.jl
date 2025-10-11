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
```@repl
using Random
A = MatrixBandwidth.random_banded_matrix(6, 1; p=1, rng=MersenneTwister(1228))
bandwidth(A)
```

Generate a ``7×7`` matrix with bandwidth ``3`` and band density ``0.3``:
```@repl
using Random
A = MatrixBandwidth.random_banded_matrix(7, 3; p=0.3, rng=MersenneTwister(0402))
bandwidth(A)
```

Generate an ``8×8`` diagonal (bandwidth ``0``) matrix with default band density (``0.5``):
```@repl
using Random
A = MatrixBandwidth.random_banded_matrix(8, 0; rng=MersenneTwister(0102))
bandwidth(A)
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
```@repl
using Graphs
g = complement(complete_multipartite_graph([3, 4, 2]))
A = Bool.(adjacency_matrix(g))
MatrixBandwidth.connected_components(A)
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
```@repl
using Graphs
g = ladder_graph(5)
A = Bool.(adjacency_matrix(g))
Int.(MatrixBandwidth.floyd_warshall_shortest_paths(A))
```

Floyd–Warshall assigns `Inf` to pairs of nodes in different connected components:
```@repl
using Graphs
g = complement(wheel_graph(8))
A = Bool.(adjacency_matrix(g))
MatrixBandwidth.floyd_warshall_shortest_paths(A)
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
```@repl
A = [4 0 9 -2; 0 0 1 0; 3 -1 5 0; 4 0 0 3]
MatrixBandwidth.is_structurally_symmetric(A)
B = [1.12 2.36 0.00; 5.99 0.0 0.0; 0.0 3.1 -7.49]
MatrixBandwidth.is_structurally_symmetric(B)
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
```@repl
A = [0 2 0 7; 0 -8 0 3; -1 9 0 0; 0 0 0 5]
MatrixBandwidth.offdiag_nz_support(A)
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
```@repl
abstract type Parent end
abstract type Child1 <: Parent end
abstract type Grandchild1 <: Child1 end
struct Grandchild2 <: Child1 end
abstract type Child2 <: Parent end
struct Child3 <: Parent end
MatrixBandwidth.find_direct_subtype(Parent, Child1)
MatrixBandwidth.find_direct_subtype(Parent, Grandchild1)
MatrixBandwidth.find_direct_subtype(Parent, Grandchild2)
MatrixBandwidth.find_direct_subtype(Parent, Child2)
MatrixBandwidth.find_direct_subtype(Parent, Child3)
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
