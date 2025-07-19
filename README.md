# MatrixBandwidth.jl

<table>
  <tr>
    <td>Metadata</td>
    <td>
      <img src="https://img.shields.io/badge/version-v0.1.0--beta-pink.svg" alt="Version">
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

## Overview

*MatrixBandwidth.jl* offers several algorithms for matrix bandwidth minimization and matrix bandwidth recognition.

The *bandwidth* of an *n*&times;*n* matrix *A* is the minimum non-negative integer *k* &isin; [0, *n* - 1] such that *A<sub>i,j</sub>* = 0 whenever |*i* - *j*| > *k*. Equivalently, *A* has bandwidth *at most* *k* if all entries above the *k*<sup>th</sup> superdiagonal and below the *k*<sup>th</sup> subdiagonal are zero, and *A* has bandwidth *at least* *k* if there exists any nonzero entry in the *k*<sup>th</sup> superdiagonal or subdiagonal.

The *matrix bandwidth minimization problem* involves finding a permutation matrix *P* such that the bandwidth of *PAP*<sup>T</sup> is minimized; this is known to be NP-complete. Several heuristic algorithms (such as Gibbs&ndash;Poole&ndash;Stockmeyer) run in polynomial time while still producing near-optimal orderings in practice, but exact methods (like Caprara&ndash;Salazar-González) are at least exponential in time complexity and thus are only feasible for relatively small matrices.

On the other hand, the *matrix bandwidth recognition problem* entails determining whether there exists a permutation matrix *P* such that the bandwidth of *PAP*<sup>T</sup> is at most some fixed non-negative integer *k* &isin; **N**&mdash;an optimal permutation that fully minimizes the bandwidth of *A* is not required. Unlike the NP-hard minimization problem, this is decidable in *O*(*n*<sup>*k*</sup>) time.

## Algorithms

The following algorithms are currently supported:

- **Minimization**
  - *Exact*
    - Caprara&ndash;Salazar-González algorithm [**under development**]
    - Del Corso&ndash;Manzini algorithm
    - Del Corso&ndash;Manzini algorithm with perimeter search
    - Saxe&ndash;Gurari&ndash;Sudborough algorithm [**under development**]
    - Brute-force search
  - *Heuristic*
    - Gibbs&ndash;Poole&ndash;Stockmeyer algorithm
    - Cuthill&ndash;McKee algorithm
    - Reverse Cuthill&ndash;McKee algorithm
  - *Metaheuristic*
    - Greedy randomized adaptive search procedure (GRASP) [**under development**]
    - Simulated annealing [**under development**]
    - Genetic algorithm [**under development**]
- **Recognition**
  - Caprara&ndash;Salazar-González algorithm [**under development**]
  - Del Corso&ndash;Manzini algorithm
  - Del Corso&ndash;Manzini algorithm with perimeter search
  - Saxe&ndash;Gurari&ndash;Sudborough algorithm [**under development**]
  - Brute-force search

(As we remain in the early stages of development, some of these may not yet be fully implemented and/or tested. Whenever an unimplemented algorithm is used, an `ERROR: TODO: Not yet implemented` is raised.)

## Installation

The only prerequisite is a working Julia installation (v1.10 or later). First, enter Pkg mode by typing `]` in the Julia REPL, then run the following command:

```julia-repl
pkg> add https://github.com/Luis-Varona/MatrixBandwidth.jl
```

When *MatrixBandwidth.jl* is finally added to the official Julia registry, you will be able to install it more easily with:

```julia-repl
pkg> add MatrixBandwidth
```

## Basic use

*MatrixBandwidth.jl* offers unified interfaces for both bandwidth minimization and bandwidth recognition via the `minimize_bandwidth` and `has_bandwidth_k_ordering` functions, respectively&mdash;the algorithm itself is specified as an argument. For example, to minimize the bandwidth of a random matrix with the reverse Cuthill&ndash;McKee algorithm, you can run the following code:

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

Similarly, to determine whether said matrix has bandwidth *at most*, say, 10 (not necessarily caring about the true minimum) via the Del Corso&ndash;Manzini algorithm, you can run:

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

If no algorithm is explicitly specified, `minimize_bandwidth` defaults to the Gibbs&ndash;Poole&ndash;Stockmeyer algorithm:

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

(We default to Gibbs&ndash;Poole&ndash;Stockmeyer because it is one of the most accurate heuristic algorithms&mdash;note how in this case, it produced a lower-bandwidth ordering than reverse Cuthill&ndash;McKee. Of course, if true optimality is required, an exact algorithm such as Caprara&ndash;Salazar-González should be used instead.)

`has_bandwidth_k_ordering` similarly defaults to Caprara&ndash;Salazar-González, which we have not yet implemented, so users should specify which of the completed algorithms they wish to use in the meantime or else face an error:

```julia-repl
julia> res_recognize_default = has_bandwidth_k_ordering(A, 10)
ERROR: TODO: Not yet implemented
[...]
```

## Documentation

The full documentation is available at [GitHub Pages](https://luis-varona.github.io/MatrixBandwidth.jl/dev/). Documentation for methods and types is also available via the Julia REPL. To learn more about the `minimize_bandwidth` function, for instance, enter help mode by typing `?`, then run the following command:

```julia-repl
help?> minimize_bandwidth
search: minimize_bandwidth MatrixBandwidth bandwidth

  minimize_bandwidth(A, solver=GibbsPooleStockmeyer()) -> MinimizationResult

  Minimize the bandwidth of A using the algorithm defined by solver.

  The bandwidth of an n×n matrix A is the minimum non-negative integer k ∈ [0,
  n - 1] such that A[i, j] = 0 whenever |i - j| > k. Equivalently, A has
  bandwidth at most k if all entries above the kᵗʰ superdiagonal and below the
  kᵗʰ subdiagonal are zero, and A has bandwidth at least k if there exists any
  nonzero entry in the kᵗʰ superdiagonal or subdiagonal.

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
  …
```

## Citing

I encourage you to cite this work if you have found any of the algorithms herein useful for your research. Starring the *MatrixBandwidth.jl* repository on GitHub is also appreciated.

The latest citation information may be found in the [CITATION.bib](https://raw.githubusercontent.com/Luis-Varona/MatrixBandwidth.jl/main/CITATION.bib) file within the repository.

## Project status

I aim to release the first stable version of *MatrixBandwidth.jl* in late July 2025. The current version is a work-in-progress, with much of the API still under development.
