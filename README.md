<img src="https://github.com/Luis-Varona/MatrixBandwidth.jl/raw/main/docs/src/assets/logo.png" alt="MatrixBandwidth.jl logo by Rebekka Jonasson" height="200"/>
<sup>MatrixBandwidth.jl logo by <a href="https://github.com/RebekkaJonasson">Rebekka Jonasson</a></sup>

# MatrixBandwidth.jl

<table>
  <tr>
    <td>Metadata</td>
    <td>
      <img src="https://img.shields.io/badge/version-v0.2.2-pink.svg" alt="Version">
      <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-A31F34.svg" alt="License: MIT"></a>
      <a href="https://github.com/JuliaDiff/BlueStyle"><img src="https://img.shields.io/badge/code%20style-blue-4495d1.svg" alt="Code Style: Blue"></a>
      <a href="https://joss.theoj.org/papers/9677739548a4f05ab398b427df1c78b9"><img src="https://joss.theoj.org/papers/9677739548a4f05ab398b427df1c78b9/status.svg" alt="JOSS: Under Review"></a>
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

MatrixBandwidth.jl offers fast algorithms for matrix bandwidth minimization and recognition. The *bandwidth* of an *n*&times;*n* matrix *A* is the minimum non-negative integer *k* &isin; {0, 1, &hellip;, *n* - 1} such that *A<sub>i,j</sub>* = 0 whenever |*i* - *j*| > *k*. Reordering the rows and columns of a matrix to reduce its bandwidth has many practical applications in engineering and scientific computing: it can improve performance when solving linear systems, approximating partial differential equations, optimizing circuit layout, and more. There are two variants of this problem: *minimization*, which involves finding a permutation matrix *P* such that the bandwidth of *PAP*<sup>T</sup> is minimized, and *recognition*, which entails determining whether there exists a permutation matrix *P* such that the bandwidth of *PAP*<sup>T</sup> is less than or equal to some fixed non-negative integer (an optimal permutation that fully minimizes the bandwidth of *A* is not required).

Many matrix bandwidth reduction algorithms exist in the literature, but implementations in the open-source ecosystem are scarce, with those that do exist primarily tackling older, less efficient algorithms. The [Boost](https://www.boost.org/) libraries in C++, the [NetworkX](https://networkx.org/) library in Python, and the MATLAB standard library all only implement the reverse Cuthill&ndash;McKee algorithm from 1971. This gap in the ecosystem not only makes it difficult for theoretical researchers to benchmark and compare new proposed algorithms but also precludes the application of the most performant modern algorithms in real-life industry settings. MatrixBandwidth.jl aims to bridge this gap by presenting a unified interface for matrix bandwidth reduction algorithms in Julia, designed with extensibility to further methods in mind.

## Algorithms

The following matrix bandwidth reduction algorithms are currently available:

- **Minimization**
  - *Exact*
    - Del Corso&ndash;Manzini
    - Del Corso&ndash;Manzini with perimeter search
    - Caprara&ndash;Salazar-González
    - Saxe&ndash;Gurari&ndash;Sudborough
    - Brute-force search
  - *Heuristic*
    - Gibbs&ndash;Poole&ndash;Stockmeyer
    - Cuthill&ndash;McKee
    - Reverse Cuthill&ndash;McKee
- **Recognition**
  - Del Corso&ndash;Manzini
  - Del Corso&ndash;Manzini with perimeter search
  - Caprara&ndash;Salazar-González
  - Saxe&ndash;Gurari&ndash;Sudborough
  - Brute-force search

Recognition algorithms determine whether any row-and-column permutation of a matrix induces bandwidth less than or equal to some fixed integer. Exact minimization algorithms always guarantee optimal orderings to minimize bandwidth, while heuristic minimization algorithms produce near-optimal solutions more quickly. Metaheuristic minimization algorithms employ iterative search frameworks to find better solutions than heuristic methods (albeit more slowly); no such algorithms are already implemented, but several are currently under development:

- Greedy randomized adaptive search procedure (GRASP)
- Particle swarm optimization with hill climbing (PSO-HC)
- Simulated annealing
- Genetic algorithm
- Ant colony optimization
- Tabu search

An index of all available algorithms by submodule (not including the unfinished metaheuristic algorithms) may also be accessed via the `MatrixBandwidth.ALGORITHMS` constant; simply run the following command in the Julia REPL:

```julia-repl
julia> MatrixBandwidth.ALGORITHMS
Dict{Symbol, Union{Dict{Symbol}, Vector}} with 2 entries:
[...]
```

If you wish to extend the interface with a new matrix bandwidth reduction algorithm, please refer to the [CONTRIBUTING.md](https://raw.githubusercontent.com/Luis-Varona/MatrixBandwidth.jl/main/CONTRIBUTING.md) file for detailed instructions.

## Installation

The only prerequisite is a working Julia installation (v1.10 or later). First, enter Pkg mode by typing `]` in the Julia REPL, then run the following command:

```julia-repl
pkg> add MatrixBandwidth
```

## Basic use

MatrixBandwidth.jl offers unified interfaces for both bandwidth minimization and bandwidth recognition via the `minimize_bandwidth` and `has_bandwidth_k_ordering` functions, respectively&mdash;the algorithm itself is specified as an argument. Comprehensive documentation for these two methods is available [here](https://luis-varona.github.io/MatrixBandwidth.jl/stable/public_api/#MatrixBandwidth.Minimization.minimize_bandwidth) and [here](https://luis-varona.github.io/MatrixBandwidth.jl/stable/public_api/#MatrixBandwidth.Recognition.has_bandwidth_k_ordering), and further details about their respective output structs can be found [here](http://luis-varona.github.io/MatrixBandwidth.jl/stable/public_api/#MatrixBandwidth.Minimization.MinimizationResult) and [here](https://luis-varona.github.io/MatrixBandwidth.jl/stable/public_api/#MatrixBandwidth.Recognition.RecognitionResult). We go over some basic examples below.

To minimize the bandwidth of a random matrix with, say, the reverse Cuthill&ndash;McKee algorithm, you can run the following code:

```julia-repl
julia> using Random, SparseArrays

julia> Random.seed!(8675309);

julia> A = sprand(60, 60, 0.01); A = A + A' # Ensure structural symmetry
60×60 SparseMatrixCSC{Float64, Int64} with 93 stored entries:
⎡⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠠⠀⢀⠀⠀⠒⠀⠀⠀⠀⠀⠀⡀⠨⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠀⠅⠀⠀⠀⠀⠀⠐⠀⠠⡀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠠⠀⠀⠀⠀⎥
⎢⠀⠂⠁⠁⢀⠐⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠁⠀⡀⠀⠀⠀⠀⠀⠀⠈⠁⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠄⡀⠀⠀⎥
⎢⠀⠂⠀⠢⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⢀⠀⠁⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠐⠀⠐⠀⠀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠁⠀⠀⎥
⎢⢠⠀⠀⠀⠀⢀⠀⠠⢀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠊⠀⠀⢠⠀⠀⠀⠀⠀⠠⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⣀⠀⠀⠁⠀⠀⠀⠄⠀⎥
⎢⡀⡈⠈⡀⠀⠀⠆⠀⠀⠀⠁⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠡⠀⠀⠔⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠁⠀⠀⠀⠠⠄⠁⎦

julia> res_minimize = minimize_bandwidth(A, Minimization.ReverseCuthillMcKee())
Results of Bandwidth Minimization Algorithm
 * Algorithm: Reverse Cuthill–McKee
 * Approach: heuristic
 * Minimum Bandwidth: 9
 * Original Bandwidth: 51
 * Matrix Size: 60×60

julia> A[res_minimize.ordering, res_minimize.ordering]
60×60 SparseMatrixCSC{Float64, Int64} with 93 stored entries:
⎡⠠⡢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠈⠠⡢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠠⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠠⡢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠢⡀⠈⠀⠂⠀⢠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⠀⠀⠀⠣⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠢⣀⡸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠪⡢⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢊⠐⢠⡀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠲⢄⡱⢀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠪⢂⡄⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠪⡢⎦
```

Similarly, to determine whether said matrix has bandwidth *at most*, say, 3 (not necessarily caring about the true minimum) via the Saxe&ndash;Gurari&ndash;Sudborough algorithm, you can run:

```julia-repl
julia> res_recognize = has_bandwidth_k_ordering(A, 3, Recognition.SaxeGurariSudborough())
Results of Bandwidth Recognition Algorithm
 * Algorithm: Saxe–Gurari–Sudborough
 * Bandwidth Threshold k: 3
 * Has Bandwidth ≤ k Ordering: true
 * Original Bandwidth: 51
 * Matrix Size: 60×60

julia> A[res_recognize.ordering, res_recognize.ordering]
60×60 SparseMatrixCSC{Float64, Int64} with 93 stored entries:
⎡⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠠⠂⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠘⠀⡠⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠐⠊⠀⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠀⠀⠠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠢⡄⠉⡢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠪⢀⠔⢂⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠰⡊⠈⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠊⠀⠠⡀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠢⣀⡸⢀⡀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠰⡀⠈⢠⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠒⣀⡸⢂⡀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠰⡀⠈⠠⡀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠢⠊⡠⎦
```

(In this case, though, it turns out that 3 *is* the true minimum bandwidth of the matrix, as can be verified by running `minimize_bandwidth` with any exact algorithm.)

If no algorithm is explicitly specified, `minimize_bandwidth` defaults to the Gibbs&ndash;Poole&ndash;Stockmeyer algorithm:

```julia-repl
julia> res_minimize_default = minimize_bandwidth(A)
Results of Bandwidth Minimization Algorithm
 * Algorithm: Gibbs–Poole–Stockmeyer
 * Approach: heuristic
 * Minimum Bandwidth: 5
 * Original Bandwidth: 51
 * Matrix Size: 60×60

julia> A[res_minimize_default.ordering, res_minimize_default.ordering]
60×60 SparseMatrixCSC{Float64, Int64} with 93 stored entries:
⎡⠠⡢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠈⠠⡢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠠⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠠⡢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠄⠁⠀⢢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⣀⣀⠘⠄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠡⠄⡡⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢⣀⠘⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⠀⡠⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠪⢂⡄⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢊⡰⢀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢤⠓⣠⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠚⠠⡢⎦
```

(We default to Gibbs&ndash;Poole&ndash;Stockmeyer because it is one of the most accurate heuristic algorithms&mdash;note how in this case, it produced a lower-bandwidth ordering than reverse Cuthill&ndash;McKee. Of course, if true optimality is required, an exact algorithm should be used instead.)

`has_bandwidth_k_ordering` similarly defaults to the Del Corso&ndash;Manzini algorithm:

```julia-repl
julia> res_recognize_default = has_bandwidth_k_ordering(A, 6)
Results of Bandwidth Recognition Algorithm
 * Algorithm: Del Corso–Manzini
 * Bandwidth Threshold k: 6
 * Has Bandwidth ≤ k Ordering: true
 * Original Bandwidth: 51
 * Matrix Size: 60×60

julia> A[res_recognize_default.ordering, res_recognize_default.ordering]
60×60 SparseMatrixCSC{Float64, Int64} with 93 stored entries:
⎡⠀⠀⠄⠑⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⢄⠁⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⢄⠀⠀⠀⠀⠠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠑⠀⡀⠀⠀⠌⠑⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠈⢆⠁⠀⠀⠀⠲⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠘⢠⡀⠐⠀⠀⡁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠄⠠⠀⠀⠐⢂⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠰⢀⡀⠈⠀⡐⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢀⠠⠀⠀⠀⣄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠀⢤⠀⠀⠀⠔⠄⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢀⠄⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠂⠀⠀⠀⠁⢀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠠⠂⢀⠑⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢄⠐⠀⠀⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎦
```

Complementing our various bandwidth minimization and recognition algorithms, MatrixBandwidth.jl exports several additional core functions, including (but not limited to) `bandwidth` and `profile` to compute the original bandwidth and profile of a matrix prior to any reordering:

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

The full documentation is available at [GitHub Pages](https://luis-varona.github.io/MatrixBandwidth.jl/). Documentation for methods and types is also available via the Julia REPL. To learn more about the `minimize_bandwidth` function, for instance, enter help mode by typing `?`, then run the following command:

```julia-repl
help?> minimize_bandwidth
search: minimize_bandwidth bandwidth MatrixBandwidth

  minimize_bandwidth(A, solver=GibbsPooleStockmeyer()) -> MinimizationResult

  Minimize the bandwidth of A using the algorithm defined by solver.

  The bandwidth of an n×n matrix A is the minimum non-negative integer k ∈
  \{0, 1, …, n - 1\} such that A[i, j] = 0 whenever |i - j| > k.

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

We encourage you to cite this work if you have found any of the algorithms herein useful for your research. Starring the MatrixBandwidth.jl repository on GitHub is also appreciated.

The latest citation information may be found in the [CITATION.bib](https://raw.githubusercontent.com/Luis-Varona/MatrixBandwidth.jl/main/CITATION.bib) file within the repository.

## Contributing

We welcome all bug reports, feature requests, and contributions! Please refer to the [CONTRIBUTING.md](https://raw.githubusercontent.com/Luis-Varona/MatrixBandwidth.jl/main/CONTRIBUTING.md) file for guidelines on how to open issues and submit pull requests.

## Project status

The latest stable release of MatrixBandwidth.jl is v0.2.2. Although several metaheuristic algorithms are still under development, the rest of the package is fully functional and covered by unit tests.

Currently, MatrixBandwidth.jl's core functions generically accept any input of the type `AbstractMatrix{<:Number}`, not behaving any differently when given sparsely stored matrices (e.g., from the [SparseArrays.jl](https://github.com/JuliaSparse/SparseArrays.jl) standard library package). Capabilities for directly handling graph inputs (aiming to reduce the matrix bandwidth of a graph's adjacency) are also not available. Given that bandwidth reduction is often applied to sparse matrices and graphs, these limitations will be addressed in a future release of the package.
