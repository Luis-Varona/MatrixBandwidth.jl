# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CapraraSalazarGonzalez <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

The *Caprara–Salazar-González minimization algorithm* is an exact method for minimizing the
bandwidth of a structurally symmetric matrix ``A``. For a fixed ``k ∈ ℕ``, the algorithm
performs a bidirectional depth-first search of all partial orderings of the rows and columns
of ``A``, adding indices one at a time to both the left and right ends. Partial orderings
are pruned not only by ensuring that adjacent pairs of currently placed indices are within
``k`` of each other but also by employing a branch-and-bound framework with lower bounds on
bandwidtth compatibility computed via integer linear programming relaxations. This search is
repeated with incrementing values of ``k`` until a bandwidth-``k`` ordering is found
[CS05], with ``k`` initialized to some lower bound on the minimum bandwidth of ``A`` up to
symmetric permutation.

Specifically, this implementation of the Caprara–Salazar-González algorithm uses the
``min(α(A), γ(A))`` lower bound from the original paper [CS05, pp. 359--60] as the initial
value of ``k``. (Further implementation details can be found in the source code for
[`bandwidth_lower_bound`](@ref).)

As noted above, the Caprara–Salazar-González algorithm requires structurally symmetric input
(that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for
``1 ≤ i, j ≤ n``).

# Supertype Hierarchy
`CapraraSalazarGonzalez` <: [`ExactSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

# Performance
[TODO: Write here]

# Examples
[TODO: Write here]

# References
- [CS05](@cite): A. Caprara and J.-J. Salazar-González. *Laying Out Sparse Graphs with
    Provably Minimum Bandwidth*. INFORMS Journal on Computing **17**, 356–73 (2005).
    https://doi.org/10.1287/ijoc.1040.0083.
"""
struct CapraraSalazarGonzalez <: ExactSolver end

# push!(ALGORITHMS[:Minimization][:Exact], CapraraSalazarGonzalez)

Base.summary(::CapraraSalazarGonzalez) = "Caprara–Salazar-González"

_requires_symmetry(::CapraraSalazarGonzalez) = true

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::CapraraSalazarGonzalez)
    error("TODO: Not yet implemented")
    return nothing
end
