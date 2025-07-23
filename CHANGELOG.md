# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Created `CHANGELOG.md` to document changes to this project.
- Clarified certain `if-else` checks in the `bandwidth` method and in a helper function for the Del Corso&ndash;Manzini `Recognition` deciders by explaining via inline comments that we cannot reduce over an empty collection.

### Changed

- Eliminated unnecessary reallocation of a boolean matrix in the `bandwidth` method by directly using `findall(!iszero, A)` instead of calling `_offdiag_nonzero_support(A)`.
- Switched from a generator comprehension in the `bandwidth` method to `Iterators.map` for more idiomatic functional programming.
- Changed the GitHub Pages references in `README.md` and the `MatrixBandwidth` module docstring from the development docs to the default site URL (which redirects to the stable documentation).
- Changed the BibTeX key in `CITATION.bib` from `Var2025` to `Var25`.
- Updated the target date for completion of core API development from mid-August 2025 to September 2025 in `README.md`.

## [0.1.0] - 2025-07-19

### Added

- Released the initial stable version and added the package to Julia's General package registry.

[unreleased]: https://github.com/Luis-Varona/MatrixBandwidth.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.0
