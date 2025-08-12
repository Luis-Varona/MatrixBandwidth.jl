# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added **References** sections to docstrings for immediate readability in the REPL and in the source code without needing to open the Documenter-generated website (#105).

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

[unreleased]: https://github.com/Luis-Varona/MatrixBandwidth.jl/compare/v0.1.3...HEAD
[0.1.3]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.3
[0.1.2]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.2
[0.1.1]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.0
