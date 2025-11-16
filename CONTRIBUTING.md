# Contributing

Welcome to MatrixBandwidth.jl! Thank you for your interest in contributing to this package. Kindly peruse these guidelines before opening an issue or submitting a pull request (PR).

## Issues

Most issues are either bug reports or feature requests; we describe the expectations for each type below. If your issue does not fit into either category, feel free to open it anyway&mdash;just be clear, concise, and of course respectful!

### Bug reports

If you encounter a bug, please include a [minimal reproducible example](https://en.wikipedia.org/wiki/Minimal_reproducible_example) containing the code necessary to reproduce the bug, along with the expected and actual behavior. Additionally, please provide information about your operating system, Julia version, and MatrixBandwidth.jl version. (If you are unsure about how to find this information, feel free to reach out anyway, and we can guide you through the process.)

### Feature requests

When suggesting a new feature, please provide a clear and concise description of your idea, along with concrete motivation for its inclusion in the package. If applicable, include references to relevant literature or existing implementations in other software. And of course, if you wish to implement the feature yourself, please indicate so in your issue, and we can get you underway to submitting a PR!

## Pull requests

Before writing a PR, we would prefer if you first reached out via an issue to check whether your contribution would fit well with the repository. Once in the process of actually writing a PR, please adhere to the following guidelines:

### Correctness

When introducing new features to the package, unit tests are of paramount importance to ensure that your code behaves as expected. Please add them to the [`test/`](test/) directory, following the existing structure and conventions: a unit test meant for a method in `src/<path>/<file>.jl` should be placed in `test/<path>/<file>.jl`.

### Documentation

Equally crucial is proper documentation. Judicious use of inline comments in complex sections is encouraged to enhance code readability, regardless of whether the method/struct is public or private. For public methods/structs, please also include docstrings matching the style of existing documentation. At bare minimum, methods should have a one-liner summary of their purpose along with descriptions of all arguments and return values, and structs should have a one-liner summary of their purpose along with descriptions of all fields. Less trivial methods/structs should have more extensive documentation&mdash;examples in particular are highly recommended.

### Style

Format your code with [JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl) in the package's environment before pushing; our [`.JuliaFormatter.toml`](.JuliaFormatter.toml) file will automatically enforce Invenia's [BlueStyle](https://github.com/invenia/BlueStyle), with a few modifications.

### Efficiency

While not always necessary, strive to write efficient code, especially for performance-critical sections. Moreover, try to avoid adding unncessary dependencies to the package; if you absolutely *must*, please justify its inclusion in your PR description.

### New algorithms

One of the most valued ways to contribute is by adding a new matrix bandwidth reduction algorithm to the interface (although bug fixes, documentation improvements, and other feature additions are always welcome). Of course, not all algorithms are suitable for inclusion, so it is best to first open an issue to discuss the possibility. While doing so, be sure to specify to which category your algorithm belongs&mdash;minimization (exact, heuristic, or metaheuristic) or recognition&mdash;and provide reference(s) to the original publication(s) describing the algorithm.

#### Adding a minimization algorithm

If your proposal is approved and you are implementing a bandwidth minimization algorithm, navigate to the appropriate directory ([`src/Minimization/Exact/`](src/Minimization/Exact/), [`src/Minimization/Heuristic/`](src/Minimization/Heuristic/), or [`src/Minimization/Metaheuristic/`](src/Minimization/Metaheuristic/)) and follow these steps:

1. Create a new file in the `solvers/` subdirectory named after your algorithm (e.g., [`cuthill_mckee.jl`](src/Minimization/Heuristic/solvers/cuthill_mckee.jl) for the Cuthill&ndash;McKee and reverse Cuthill&ndash;McKee algorithms).
2. Define a new subtype of `ExactSolver`, `HeuristicSolver`, or `MetaheuristicSolver` (as appropriate) to represent your algorithm. For instance, the reverse Cuthill&ndash;McKee algorithm uses the `ReverseCuthillMcKee` struct in [`solvers/cuthill_mckee.jl`](src/Minimization/Heuristic/solvers/cuthill_mckee.jl).
3. Add the algorithm to the appropriate entry of the `MatrixBandwidth.ALGORITHMS` constant, like so: `push!(MatrixBandwidth.ALGORITHMS[:Minimization][:Heuristic], ReverseCuthillMcKee)`.
4. Overload `Base.summary` to return the name of your algorithm and `MatrixBandwidth._requires_structural_symmetry` to indicate whether your algorithm requires a structurally symmetric matrix. For example: `Base.summary(::ReverseCuthillMcKee) = "Reverse Cuthill–McKee"` and `MatrixBandwidth._requires_structural_symmetry(::ReverseCuthillMcKee) = true`.
5. Implement a corresponding `Minimization._minimize_bandwidth_impl(::AbstractMatrix{Bool}, ::NewSolverType)` method with the core algorithm logic. Do *not* attempt to directly overload the `minimize_bandwidth` function, as this contains common preprocessing, postprocessing, and input validation logic shared by all minimization algorithms.
6. Export your new solver struct from both [`src/Minimization/Minimization.jl`](src/Minimization/Minimization.jl) and the appropriate submodule file (either [`Minimization/Exact/Exact.jl`](src/Minimization/Exact/Exact.jl), [`Minimization/Heuristic/Heuristic.jl`](src/Minimization/Heuristic/Heuristic.jl), or [`Minimization/Metaheuristic/Metaheuristic.jl`](src/Minimization/Metaheuristic/Metaheuristic.jl)).
7. Properly document your solver struct and add unit tests in the corresponding test file, adhering to the guidelines outlined above.

For a good reference implementation, you may look at [`src/Minimization/Heuristic/solvers/cuthill_mckee.jl`](src/Minimization/Heuristic/solvers/cuthill_mckee.jl) and its associated test suite in [`test/Minimization/Heuristic/solvers/cuthill_mckee.jl`](test/Minimization/Heuristic/solvers/cuthill_mckee.jl). Also refer to [`src/Minimization/Minimization.jl`](src/Minimization/Minimization.jl) and [`src/Minimization/Heuristic/Heuristic.jl`](src/Minimization/Heuristic/Heuristic.jl) for guidance on exporting your new solver struct.

#### Adding a recognition algorithm

The process for adding a bandwidth recognition algorithm, once your proposal is approved, is very similar. First, navigate to the [`src/Recognition/`](src/Recognition/) directory, then follow these steps:

1. Create a new file in the [`deciders/`](src/Recognition/deciders/) subdirectory named after your algorithm (e.g., [`del_corso_manzini.jl`](src/Recognition/deciders/del_corso_manzini.jl) for the Del Corso&ndash;Manzini and Del Corso&ndash;Manzini with perimeter search algorithms).
2. Define a new subtype of `AbstractDecider` to represent your algorithm. For instance, the Del Corso&ndash;Manzini algorithm uses the `DelCorsoManzini` struct in [`deciders/del_corso_manzini.jl`](src/Recognition/deciders/del_corso_manzini.jl).
3. Add the algorithm to the `MatrixBandwidth.ALGORITHMS[:Recognition]` constant, like so: `push!(MatrixBandwidth.ALGORITHMS[:Recognition], DelCorsoManzini)`.
4. Overload `Base.summary` to return the name of your algorithm and `MatrixBandwidth._requires_structural_symmetry` to indicate whether your algorithm requires a structurally symmetric matrix. For example: `Base.summary(::DelCorsoManzini) = "Del Corso–Manzini"` and `MatrixBandwidth._requires_structural_symmetry(::DelCorsoManzini) = true`.
5. Implement a corresponding `Recognition._has_bandwidth_k_ordering_impl(::AbstractMatrix{Bool}, ::Integer, ::NewDeciderType)` method with the core algorithm logic. Do *not* attempt to directly overload the `has_bandwidth_k_ordering` function, as this contains common preprocessing, postprocessing, and input validation logic shared by all recognition algorithms.
6. Export your new decider struct from [`src/Recognition/Recognition.jl`](src/Recognition/Recognition.jl).
7. Properly document your decider struct and add unit tests in the corresponding test file, adhering to the guidelines outlined above.

For a good reference implementation, you may look at [`src/Recognition/deciders/del_corso_manzini.jl`](src/Recognition/deciders/del_corso_manzini.jl) and its associated test suite in [`test/Recognition/deciders/del_corso_manzini.jl`](test/Recognition/deciders/del_corso_manzini.jl). Also refer to [`src/Recognition/Recognition.jl`](src/Recognition/Recognition.jl) for guidance on exporting your new decider struct.

### Citing

When adding substantial new features, it is important to cite relevant literature in the documentation. Add a `# References` section in the docstring of the relevant method/struct, and include citations in [BibTeX](https://en.wikipedia.org/wiki/BibTeX) format in [`docs/src/refs.bib`](docs/src/refs.bib). For example, the `profile` function in [`src/core.jl`](src/core.jl) includes the citation

```julia
"""
# References
- [Maf14](@cite): L. O. Mafteiu-Scai. *The Bandwidths of a Matrix. A Survey of Algorithms*.
    Annals of West University of Timisoara - Mathematics and Computer Science **52**,
    183–223 (2014). https://doi.org/10.2478/awutm-2014-0019.
"""
```

with corresponding BibTeX entry

```bibtex
@article{Maf14,
  author = {Mafteiu-Scai, Liviu Octavian},
  journal = {Annals of West University of Timisoara - Mathematics and Computer Science},
  number = {2},
  title = {The Bandwidths of a Matrix. A Survey of Algorithms},
  volume = {52},
  year = {2014},
  pages = {183--223},
  url = {https://doi.org/10.2478/awutm-2014-0019}
}
```
