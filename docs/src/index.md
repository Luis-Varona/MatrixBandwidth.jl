```@meta
CurrentModule = MatrixBandwidth
```

```@raw html
<table>
  <tr>
    <td>Metadata</td>
    <td>
      <img src="https://img.shields.io/badge/version-v0.1.3-pink.svg" alt="Version">
      <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-A31F34.svg" alt="License: MIT"></a>
      <a href="https://github.com/JuliaDiff/BlueStyle"><img src="https://img.shields.io/badge/code%20style-blue-4495d1.svg" alt="Code Style: Blue"></a>
    </td>
  </tr>
  <tr>
    <td>Documentation</td>
    <td>
      <a href="https://luis-varona.github.io/MatrixBandwidth.jl/stable/"><img src="https://img.shields.io/badge/docs-stable-darkgreen.svg" alt="Documentation of latest stable version"></a>
      <a href="https://luis-varona.github.io/MatrixBandwidth.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-rebeccapurple.svg" alt="Documentation of dev version"></a>
    </td>
  </tr>
  <tr>
    <td>Continuous integration</td>
    <td>
      <a href="https://github.com/Luis-Varona/MatrixBandwidth.jl/actions?query=workflow%3ACI+branch%3Amain"><img src="https://github.com/Luis-Varona/MatrixBandwidth.jl/actions/workflows/CI.yml/badge.svg" alt="GitHub Workflow Status"></a>
    </td>
  </tr>
  <tr>
    <td>Code coverage</td>
    <td>
      <a href="https://codecov.io/gh/Luis-Varona/MatrixBandwidth.jl"><img src="https://codecov.io/gh/Luis-Varona/MatrixBandwidth.jl/branch/main/graph/badge.svg" alt="Test coverage from codecov"></a>
    </td>
    </tr>
    <tr>
      <td>Static analysis with</td>
      <td>
        <a href="https://github.com/JuliaTesting/Aqua.jl"><img src="https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg" alt="Aqua QA"></a>
        <a href="https://github.com/aviatesk/JET.jl"><img src="https://img.shields.io/badge/%E2%9C%88%20tested%20with-JET.jl%EF%B8%8F-9cf.svg" alt="JET static analysis"></a>
      </td>
    </tr>
</table>
```

## Overview

*MatrixBandwidth.jl* offers fast algorithms for matrix bandwidth minimization and matrix bandwidth recognition. Reordering the rows and columns of a matrix to reduce its bandwidth has many practical applications in engineering and scientific computing. It is a common preprocessing step used to improve performance when solving linear systems, approximating partial differential equations, optimizing circuit layout, and more.

Recall that the *bandwidth* of an ``n \times n`` matrix ``A`` is the minimum non-negative integer ``k \in \{0, 1, \ldots, n - 1\}`` such that ``A_{i,j} = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most* ``k`` if all entries above the ``k^\text{th}`` superdiagonal and below the ``k^\text{th}`` subdiagonal are zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the ``k^\text{th}`` superdiagonal or subdiagonal.

The *matrix bandwidth minimization problem* involves finding a permutation matrix ``P`` such that the bandwidth of ``PAP^\mathsf{T}`` is minimized; this is known to be NP-complete. Several heuristic algorithms (such as Gibbs–Poole–Stockmeyer) run in polynomial time while still producing near-optimal orderings in practice, but exact methods (like Caprara–Salazar-González) are at least exponential in time complexity and thus are only feasible for relatively small matrices.

On the other hand, the *matrix bandwidth recognition problem* entails determining whether there exists a permutation matrix ``P`` such that the bandwidth of ``PAP^\mathsf{T}`` is at most some fixed non-negative integer ``k \in \mathbb{N}``—an optimal permutation that fully minimizes the bandwidth of ``A`` is not required. Unlike the NP-hard minimization problem, this is decidable in ``O(n^k)`` time.

## Algorithms

The following algorithms are currently supported:

- **Minimization**
  - *Exact*
    - Caprara–Salazar-González algorithm [**under development**]
    - Del Corso–Manzini algorithm
    - Del Corso–Manzini algorithm with perimeter search
    - Saxe–Gurari–Sudborough algorithm [**under development**]
    - Brute-force search
  - *Heuristic*
    - Gibbs–Poole–Stockmeyer algorithm
    - Cuthill–McKee algorithm
    - Reverse Cuthill–McKee algorithm
  - *Metaheuristic*
    - Greedy randomized adaptive search procedure (GRASP) [**under development**]
    - Simulated annealing [**under development**]
    - Genetic algorithm [**under development**]
- **Recognition**
  - Caprara–Salazar-González algorithm [**under development**]
  - Del Corso–Manzini algorithm
  - Del Corso–Manzini algorithm with perimeter search
  - Saxe–Gurari–Sudborough algorithm [**under development**]
  - Brute-force search

(Although the API is already stable with the bulk of the library already functional and tested, a few algorithms remain under development. Whenever such an algorithm is used, the error `ERROR: TODO: Not yet implemented` is raised.)

An index of all available algorithms by submodule may also be accessed via the `MatrixBandwidth.ALGORITHMS` constant; simply run the following command in the Julia REPL:

```julia-repl
julia> MatrixBandwidth.ALGORITHMS
Dict{Symbol, Union{Dict{Symbol}, Vector}} with 2 entries:
[...]
```

## Installation

The only prerequisite is a working Julia installation (v1.10 or later). First, enter Pkg mode by typing `]` in the Julia REPL, then run the following command:

```julia-repl
pkg> add MatrixBandwidth
```

## Basic use

*MatrixBandwidth.jl* offers unified interfaces for both bandwidth minimization and bandwidth recognition via the `minimize_bandwidth` and `has_bandwidth_k_ordering` functions, respectively—the algorithm itself is specified as an argument. For example, to minimize the bandwidth of a random matrix with the reverse Cuthill–McKee algorithm, you can run the following code:

```julia-repl
julia> using Random, SparseArrays

julia> Random.seed!(8675309);

julia> A = sprand(40, 40, 0.02); A = A + A' # Ensure structural symmetry
40×40 SparseMatrixCSC{Float64, Int64} with 82 stored entries:
⎡⠀⠀⠂⠘⠀⠀⠐⠀⠀⠀⠀⠀⠆⠀⠀⠀⠂⠄⠈⠀⎤
⎢⣈⠀⠤⠃⠀⠀⠀⠈⠀⠀⠀⠀⠂⠂⠀⠀⠀⠂⠀⠄⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠚⠀⠀⠀⠀⠀⠂⠁⠀⠀⠀⠀⎥
⎢⠐⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠢⠀⠈⠀⠀⠂⡄⎥
⎢⠀⠀⠀⠀⠚⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠠⠀⠒⎥
⎢⠀⠀⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠂⠨⠂⠀⠀⠁⠀⎥
⎢⠈⠁⠨⠀⠀⠀⠠⡀⠀⠀⠠⠀⠀⠀⠀⠂⠠⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠌⠀⡀⠀⠀⠀⠢⠂⠠⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠈⠄⠠⠀⠀⠀⠀⠀⠀⡐⠀⠀⠀⠂⠀⠀⢀⡰⠀⠀⎥
⎣⠂⠀⠀⠄⠀⠀⠈⠤⢠⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⎦

julia> res_minimize = minimize_bandwidth(A, Minimization.ReverseCuthillMcKee())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Reverse Cuthill–McKee
 * Approach: heuristic
 * Minimum Bandwidth: 9
 * Original Bandwidth: 37
 * Matrix Size: 40×40

julia> A[res_minimize.ordering, res_minimize.ordering]
40×40 SparseMatrixCSC{Float64, Int64} with 82 stored entries:
⎡⠪⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠀⢠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠒⡀⠈⠆⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠈⠡⠄⠁⠁⠠⠂⢄⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠁⡀⠀⠀⠠⠀⠘⢄⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠈⢄⠀⠂⢀⠐⠀⠠⠁⢆⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠒⢄⠀⡀⠀⠀⠀⠌⢢⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠡⢄⡀⠄⠀⠀⠘⡄⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠒⠒⠤⢄⡱⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠪⡢⎦
```

Similarly, to determine whether said matrix has bandwidth *at most*, say, 10 (not necessarily caring about the true minimum) via the Del Corso–Manzini algorithm, you can run:

```julia-repl
julia> res_recognize = has_bandwidth_k_ordering(A, 10, Recognition.DelCorsoManzini())
Results of Bandwidth Recognition Algorithm
 * Algorithm: Del Corso–Manzini
 * Bandwidth Threshold k: 10
 * Has Bandwidth ≤ k Ordering: true
 * Original Bandwidth: 37
 * Matrix Size: 40×40

julia> A[res_recognize.ordering, res_recognize.ordering]
40×40 SparseMatrixCSC{Float64, Int64} with 82 stored entries:
⎡⠊⠀⠀⢀⡈⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⢀⠀⠀⢀⠀⠀⡐⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⢆⠈⠀⠐⢀⠐⠀⠅⡀⠤⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠑⢀⠠⠄⠄⠊⠀⠈⠀⠀⠐⢀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠑⠀⡌⠂⠀⢀⠐⠀⠀⠀⠔⡀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠑⢀⠀⠀⠀⠀⠀⠀⠀⠐⠲⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠐⢀⠄⠀⠀⠀⠀⠀⠀⠁⢁⠄⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢰⡀⠀⠀⠀⠀⠀⠑⢀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠅⢀⢄⠀⠀⠀⠀⢀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠐⠀⢀⠀⠀⎦
```

If no algorithm is explicitly specified, `minimize_bandwidth` defaults to the Gibbs–Poole–Stockmeyer algorithm:

```julia-repl
julia> res_minimize_default = minimize_bandwidth(A)
Results of Bandwidth Minimization Algorithm
 * Algorithm: Gibbs–Poole–Stockmeyer
 * Approach: heuristic
 * Minimum Bandwidth: 6
 * Original Bandwidth: 37
 * Matrix Size: 40×40

julia> A[res_minimize_default.ordering, res_minimize_default.ordering]
40×40 SparseMatrixCSC{Float64, Int64} with 82 stored entries:
⎡⠪⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⡠⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠐⠀⠀⠑⠄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠑⠄⠀⠀⠌⢢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠈⠢⣁⠀⠀⠨⠆⢀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠈⠢⠆⠄⡡⠚⠄⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠚⠄⡀⠈⠦⣀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢣⠀⠀⠃⡀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠠⠎⡡⠢⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠂⢠⠒⎦
```

(We default to Gibbs–Poole–Stockmeyer because it is one of the most accurate heuristic algorithms—note how in this case, it produced a lower-bandwidth ordering than reverse Cuthill–McKee. Of course, if true optimality is required, an exact algorithm such as Caprara–Salazar-González should be used instead.)

`has_bandwidth_k_ordering` similarly defaults to Caprara–Salazar-González, which we have not yet implemented, so users should specify which of the completed algorithms they wish to use in the meantime or else face an error:

```julia-repl
julia> res_recognize_default = has_bandwidth_k_ordering(A, 10)
ERROR: TODO: Not yet implemented
[...]
```

Complementing our various bandwidth minimization and recognition algorithms, *MatrixBandwidth.jl* exports several additional core functions, including (but not limited to) `bandwidth` and `profile` to compute the original bandwidth and profile of a matrix:

```julia-repl
julia> using Random, SparseArrays

julia> Random.seed!(1234);

julia> A = sprand(50, 50, 0.02)
50×50 SparseMatrixCSC{Float64, Int64} with 49 stored entries:
⎡⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠂⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⡀⠀⡀⠀⠀⠀⠄⠀⠀⠀⠄⠀⎥
⎢⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠄⠂⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠈⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⢀⠀⠀⠀⠂⠀⠀⠀⠁⎥
⎢⡀⡀⢀⠄⠀⠁⠄⢀⠀⠀⠀⠀⠀⢀⠀⠀⠠⠀⠀⠀⠀⣀⠀⠀⠀⎥
⎢⠀⢀⠀⠀⠀⠊⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⎥
⎢⠈⠀⢀⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠨⠐⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎦

julia> bandwidth(A) # Bandwidth prior to any reordering of rows and columns
38

julia> profile(A) # Profile prior to any reordering of rows and columns
703
```

(Closely related to bandwidth, the *column profile* of a matrix is the sum of the distances from each diagonal entry to the farthest nonzero entry in that column, whereas the *row profile* is the sum of the distances from each diagonal entry to the farthest nonzero entry in that row. `profile(A)` computes the column profile of `A` by default, but it can also be used to compute the row profile.)

## Documentation

The full documentation is available at [GitHub Pages](https://luis-varona.github.io/MatrixBandwidth.jl/). Documentation for methods and types is also available via the Julia REPL—for instance, to learn more about the `minimize_bandwidth` function, enter help mode by typing `?`, then run the following command:

```julia-repl
help?> minimize_bandwidth
search: minimize_bandwidth bandwidth MatrixBandwidth

  minimize_bandwidth(A, solver=GibbsPooleStockmeyer()) -> MinimizationResult

  Minimize the bandwidth of A using the algorithm defined by solver.

  The bandwidth of an n×n matrix A is the minimum non-negative integer k ∈
  \{0, 1, …, n - 1\} such that A[i, j] = 0 whenever |i - j| > k. Equivalently,
  A has bandwidth at most k if all entries above the kᵗʰ superdiagonal and
  below the kᵗʰ subdiagonal are zero, and A has bandwidth at least k if there
  exists any nonzero entry in the kᵗʰ superdiagonal or subdiagonal.

  This function computes a (near-)optimal ordering π of the rows and columns
  of A so that the bandwidth of PAPᵀ is minimized, where P is the permutation
  matrix corresponding to π. This is known to be an NP-complete problem;
  however, several heuristic algorithms such as Gibbs–Poole–Stockmeyer run in
  polynomial time while still still producing near-optimal orderings in
  practice. Exact methods like Caprara–Salazar-González are also available,
  but they are at least exponential in time complexity and thus only feasible
  for relatively small matrices.

  Arguments
  ≡≡≡≡≡≡≡≡≡
  [...]
```

## Citing

I encourage you to cite this work if you have found any of the algorithms herein useful for your research. Starring the *MatrixBandwidth.jl* repository on GitHub is also appreciated.

The latest citation information may be found in the [CITATION.bib](https://raw.githubusercontent.com/Luis-Varona/MatrixBandwidth.jl/main/CITATION.bib) file within the repository.

## Project status

The latest stable release of *MatrixBandwidth.jl* is v0.1.3. Although several algorithms are still under development, the bulk of the library is already functional and tested. I aim to complete development (including documentation and tests) of the remaining algorithms and other utility features by October 2025.

## Index

```@index
```
