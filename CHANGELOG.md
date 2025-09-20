# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2025-09-20

### Added

- Declared additional unfinished metaheuristic minimization algorithms (`PSOHC`, `AntColony`, and `TabuSearch`) in `MatrixBandwidth.Minimization.Metaheuristic` (#155).
- Added doctest examples for the `minimize_bandwidth` and `has_bandwidth_k_ordering` functions (#154).
- Implemented the Caprara&ndash;Salazar-Gonz&aacute;lez algorithm (both the minimization solver and the recognition decider), concurrently adding *JuMP.jl* and *HiGHS.jl* as notable dependencies for integer linear programming (#137, #149, #152, #156).

### Changed

- Switched the default recognition decider from Caprara&ndash;Salazar-Gonz&aacute;lez to Del Corso&ndash;Manzini (#149).
- Added instructions to the documentation for extending the interface with new algorithms via `_minimize_bandwidth_impl` and `_has_bandwidth_k_ordering_impl` (#145).
- Renamed `_bool_minimal_band_ordering` and `_bool_bandwidth_k_ordering` (the functions wherein new algorithm logic is defined) to `_minimize_bandwidth_impl` and `_has_bandwidth_k_ordering_impl`, respectively (#145).
- Split up the `@index` and `@autodocs` blocks in the Documenter-generated website by submodule (#144).
- Modified `minimize_bandwidth` to skip the algorithm call and simply use the original ordering in all cases when `bandwidth(A) == bandwidth_lower_bound(A)`, not just when `bandwidth(A) == 0` (#143).
- Removed the export of `random_banded_matrix` from `MatrixBandwidth` (#140).
- Exported `AbstractAlgorithm` and `AbstractResult` from `MatrixBandwidth`, and exported their various subtypes from the appropriate submodules (#140).
- Reformatted the export blocks in each submodule to be more in cohesive, similarly to how *Graphs.jl* does it (#140).
- Removed the leading underscore from many (but not all) top-level utility functions (also renaming `_offdiag_nonzero_support` to `offdiag_nz_support`), adding docstrings to all of these (#140).
- Refactored the import structure in submodules to do, e.g., `using MatrixBandwidth` in favor of a bunch of `import ..A, ..B, ..C` statements, similarly to how *Graphs.jl* does it (#140).

### Removed

- Removed the empty test suites for the core functions of `MatrixBandwidth.Minimization` and `MatrixBandwidth.Recognition`, as `minimize_bandwidth` and `has_bandwidth_k_ordering` are already tested via the algorithm-specific test suites (#154).

## [0.1.4] - 2025-09-03

### Added

- Updated the module docstring and `README.md` to include a paragraph on the motivation and goals of this package (#131).
- Implemented the Saxe&ndash;Gurari&ndash;Sudborough algorithm (both the minimization solver and the recognition decider) (#123, #126).
- Added `bi_criteria_node_finder` (an improvement upon `pseudo_peripheral_node_finder`) as a new node finder for the heuristic solvers (#112).
- Finished unit tests for all root-level utility functions (#108, #109).
- Added **References** sections to docstrings for immediate readability in the REPL and in the source code without needing to open the Documenter-generated website (#105).

### Changed

- Replaced all deprecated `DataStructures.enqueue!` and `DataStructures.dequeue!` calls with `Base.push!` and `Base.popfirst!`, respectively (#129).
- Fleshed out and fixed a few typos in the documentation for the Del Corso&ndash;Manzini algorithms (#126).
- Renamed the `_requires_symmetry` internal function (used for input validation) to `_requires_structural_symmetry` (#123).
- Fleshed out documentation (particularly inline comments) for the Gibbs&ndash;Poole&ndash;Stockmeyer source code (#116).
- Improved unit tests for the heuristic solvers with more edge cases and scenarios (including the use of custom node finders) (#112).
- Changed the `DEFAULT_NODE_FINDER` constant for the heuristic solvers from `pseudo_peripheral_node_finder` to `bi_criteria_node_finder` (#112).
- Moved the `_connected_components` function from `MatrixBandwidth.Minimization.Heuristic` to the root `MatrixBandwidth` module (specifically `src/utils.jl`) for universal access (#109).

### Fixed

- Tightened the `bandwidth_lower_bound` function for disconnected graphs by taking the maximum of the bounds computed for each connected component (#124).
- Fixed some test names in the Del Corso&ndash;Manzini recognition algorithm test suite ("Bandwidth < k" was meant to be "Bandwidth > k", and "Bandwidth ≥ k" was meant to be "Bandwidth ≤ k") (#123).

## [0.1.3] - 2025-08-05

### Changed

- Bumped compat for *DataStructures.jl* from `0.18.15` to `0.18.15 - 0.19` (#100).

## [0.1.2] - 2025-07-31

### Added

- Started using *PrecompileTools.jl* to compile all solvers/deciders during package startup, reducing delay on first usage (formerly reached up to ~3 seconds for some algorithms) (#93).
- Created the `MatrixBandwidth.ALGORITHMS` constant to index all available algorithms by submodule (#93).
- Added "Supertype Hierarchy" sections to the docstrings of all subtypes (both concrete and abstract) (#90).
- Added a section in the `MatrixBandwidth` module docstring covering practical applications of matrix bandwidth reduction in engineering and scientific computing (#86).

### Changed

- Changed some user-facing parameters typed as `Int` to the more generic `Integer` (#90).

## [0.1.1] - 2025-07-26

### Added

- Added the `profile` function to compute the original profile of a matrix prior to any reordering (#78).
- Added the `homepage` field to `Project.toml` to reference the GitHub Pages documentation (#76).
- Created `CHANGELOG.md` to document changes to this project (#72).
- Clarified certain `if-else` checks in the `bandwidth` method and in a helper function for the Del Corso&ndash;Manzini `Recognition` deciders by explaining via inline comments that we cannot reduce over an empty collection (#71).

### Changed

- Added PR numbers to changelog entries for better traceability (#73).
- Eliminated unnecessary reallocation of a boolean matrix in the `bandwidth` method by directly using `findall(!iszero, A)` instead of calling `_offdiag_nonzero_support(A)` (#71).
- Switched from a generator comprehension in the `bandwidth` method to `Iterators.map` (more idiomatic) (#71).

### Fixed

- Changed some blocks enclosed in double backticks to be enclosed in single backticks instead (meant to be rendered as code blocks, not mathematical expressions) (#78).
- Fixed the rendering of the `dcm_ps_optimal_depth` docstring (#78).
- Updated the compatibility requirements in `test/Project.toml` to allow only a finite number of breaking releases of *Aqua.jl* and *JET.jl* (#74).

## [0.1.0] - 2025-07-19

### Added

- Released the initial stable version of the package.

[unreleased]: https://github.com/Luis-Varona/MatrixBandwidth.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.2.0
[0.1.4]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.4
[0.1.3]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.3
[0.1.2]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.2
[0.1.1]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.0
