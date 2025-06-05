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
      <a href="https://Luis-Varona.github.io/MatrixBandwidth.jl/stable/"><img src="https://img.shields.io/badge/docs-stable-darkgreen.svg" alt="Documentation of latest stable version"></a>
      <a href="https://Luis-Varona.github.io/MatrixBandwidth.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-rebeccapurple.svg" alt="Documentation of dev version"></a>
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
      <a href="https://codecov.io/gh/Luis-Varona/MatrixBandwidth.jl"><img src="https://img.shields.io/codecov/c/gh/Luis-Varona/MatrixBandwidth.jl.svg?label=codecov" alt="Test coverage from codecov"></a>
    </td>
  </tr>
</table>

## Overview

*MatrixBandwidth.jl* offers several exact, heuristic, and metaheuristic algorithms for matrix bandwidth minimization.

## Algorithms

- Exact
  - Minimum bandwidth by iterative deepening (MB-ID)
  - Minimum bandwidth by perimeter search (MB-PS)
- Heuristic
  - Cuthill&ndash;McKee algorithm
  - Reverse Cuthill&ndash;McKee algorithm
- Metaheuristic
  - Simulated annealing
  - Genetic algorithm
  - Greedy randomized adaptive search procedure (GRASP)

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

I aim to release the first stable version of *MatrixBandwidth.jl* in mid-June 2025. The current version is a work-in-progress, with much of the API still under development.
