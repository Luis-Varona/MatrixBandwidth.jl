# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added "Supertype Hierarchy" sections to the docstrings of all subtypes (both concrete and abstract) (#90).
- Added a section in `README.md` and the `MatrixBandwidth` module docstring covering practical applications of matrix bandwidth reduction in engineering and scientific computing (#86).

### Changed

- Changed some user-facing parameters typed as `Int` to the more generic `Integer` (#90).
- Changed the phrase `… integer k ∈ [0, n - 1] …` to `… integer k ∈ {0, 1, …, n - 1} …` every time matrix bandwidth is defined in the documentation (#84).

## [0.1.1] - 2025-07-26

### Added

- Mentioned the core package exports (including the new `profile` function) in the `MatrixBandwidth` module docstring and `README.md` (#78).
- Added the `profile` function to compute the original profile of a matrix prior to any reordering (#78).
- Added the `homepage` field to `Project.toml` to reference the GitHub Pages documentation (#76).
- Created `CHANGELOG.md` to document changes to this project (#72).
- Clarified certain `if-else` checks in the `bandwidth` method and in a helper function for the Del Corso&ndash;Manzini `Recognition` deciders by explaining via inline comments that we cannot reduce over an empty collection (#71).

### Changed

- Referenced the new `profile` function in the `CuthillMcKee`, `ReverseCuthillMcKee`, and `GibbsPooleStockmeyer` docstrings when matrix profile is mentioned (#78).
- Reformatted some entries in `docs/src/refs.bib` (some author names did not exactly match the papers) (#78).
- Changed "*MatrixBandwidth.jl* offers several algorithms&hellip;" to "*MatrixBandwidth.jl* offers fast algorithms&hellip;" in `README.md`. Similarly changed "Luis-Varona/MatrixBandwidth.jl: Algorithms&hellip;" to "Luis-Varona/MatrixBandwidth.jl: Fast algorithms&hellip;" in `CITATION.bib` (#74).
- Added PR numbers to changelog entries for better traceability (#73).
- Eliminated unnecessary reallocation of a boolean matrix in the `bandwidth` method by directly using `findall(!iszero, A)` instead of calling `_offdiag_nonzero_support(A)` (#71).
- Switched from a generator comprehension in the `bandwidth` method to `Iterators.map` (more idiomatic) (#71).
- Changed the GitHub Pages references in `README.md` and the `MatrixBandwidth` module docstring from the development docs to the default site URL (which redirects to the stable documentation) (#71).
- Changed the BibTeX key in `CITATION.bib` from `Var2025` to `Var25` (#71).
- Updated the target date for completion of core API development from mid-August 2025 to September 2025 in `README.md` (#71).

### Fixed

- Changed some blocks enclosed in double backticks to be enclosed in single backticks instead (meant to be rendered as code blocks, not mathematical expressions) (#78).
- Fixed the rendering of the `dcm_ps_optimal_depth` docstring (#78).
- Fixed the copyright preface in `docs/make.jl` (#75).
- Updated the compatibility requirements in `test/Project.toml` to allow only a finite number of breaking releases of `Aqua` and `JET` (#74).

## [0.1.0] - 2025-07-19

### Added

- Released the initial stable version and added the package to Julia's General package registry.

[unreleased]: https://github.com/Luis-Varona/MatrixBandwidth.jl/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.0
