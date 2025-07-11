# MatrixBandwidth.jl

<table>
  <tr>
    <td>Metadata</td>
    <td>
      <img src="https://img.shields.io/badge/version-v0.1.0--dev-pink.svg" alt="Version">
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

The *matrix bandwidth minimization problem* involves finding a permutation matrix *P* such that the bandwidth of *PAP*<sup>T</sup> is minimized; this is known to be NP-complete. Several heuristic algorithms (such as Gibbs&ndash;Poole&ndash;Stockmeyer) run in polynomial time while still producing near-optimal orderings in practice, but exact methods (like Caprara&ndash;Salazar-González) are exponential in time complexity and thus are only feasible for relatively small matrices.

On the other hand, the *matrix bandwidth recognition problem* entails determining whether there exists a permutation matrix *P* such that the bandwidth of *PAP*<sup>T</sup> is at most some fixed non-negative integer *k* &isin; **N**&mdash;an optimal permutation that fully minimizes the bandwidth of *A* is not required. Unlike the NP-hard minimization problem, this is decidable in *O*(*n*<sup>*k*</sup>) time.

## Algorithms

The following algorithms are currently supported:

- **Minimization**
  - *Exact*
    - Caprara&ndash;Salazar-González algorithm
    - Del Corso&ndash;Manzini algorithm
    - Del Corso&ndash;Manzini algorithm with perimeter search
    - Saxe&ndash;Gurari&ndash;Sudborough algorithm
  - *Heuristic*
    - Gibbs&ndash;Poole&ndash;Stockmeyer algorithm
    - Cuthill&ndash;McKee algorithm
    - Reverse Cuthill&ndash;McKee algorithm
  - *Metaheuristic*
    - Greedy randomized adaptive search procedure (GRASP)
    - Simulated annealing
    - Genetic algorithm
- **Recognition**
  - Caprara&ndash;Salazar-González algorithm
  - Del Corso&ndash;Manzini algorithm
  - Del Corso&ndash;Manzini algorithm with perimeter search
  - Saxe&ndash;Gurari&ndash;Sudborough algorithm

(As we remain in the early stages of development, some of these may not yet be fully implemented and/or tested.)

## Installation

The only prerequisite is a working Julia installation (v1.10 or later). First, enter Pkg mode by typing `]` in the Julia REPL, then run the following command:

```julia-repl
pkg> add https://github.com/Luis-Varona/MatrixBandwidth.jl
```

When *MatrixBandwidth.jl* is finally added to the official Julia registry, you will be able to install it more easily with:

```julia-repl
pkg> add MatrixBandwidth
```

## Citing

I encourage you to cite this work if you have found any of the algorithms herein useful for your research. Starring the *MatrixBandwidth.jl* repository on GitHub is also appreciated.

The latest citation information may be found in the [CITATION.bib](https://raw.githubusercontent.com/Luis-Varona/MatrixBandwidth.jl/main/CITATION.bib) file within the repository.

## Project status

I aim to release the first stable version of *MatrixBandwidth.jl* in late July 2025. The current version is a work-in-progress, with much of the API still under development.
