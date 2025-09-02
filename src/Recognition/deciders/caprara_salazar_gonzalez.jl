# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CapraraSalazarGonzalez <: AbstractDecider <: AbstractAlgorithm

The *Caprara–Salazar-González recognition algorithm* is a method for determining, given some
fixed ``k ∈ ℕ``, whether a structurally symmetric matrix ``A`` has a bandwidth at most ``k``
up to symmetric permutation. The algorithm performs a bidirectional depth-first search of
all partial orderings of the rows and columns of ``A``, adding indices one at a time to both
the left and right ends. Partial orderings are pruned not only by ensuring that adjacent
pairs of currently placed indices are within ``k`` of each other but also by employing a
branch-and-bound framework with lower bounds on bandwidtth compatibility computed via
integer linear programming relaxations. This search is repeated with incrementing values of
``k`` until a bandwidth-``k`` ordering is found [CS05], with ``k`` initialized to some lower
bound on the minimum bandwidth of ``A`` up to symmetric permutation.

As noted above, the Caprara–Salazar-González algorithm requires structurally symmetric input
(that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for
``1 ≤ i, j ≤ n``).

# Supertype Hierarchy
`CapraraSalazarGonzalez` <: [`AbstractDecider`](@ref) <: [`AbstractAlgorithm`](@ref)

# Performance
[TODO: Write here]

# Examples
[TODO: Write here]

# Notes
This algorithm is not the main one described in the original paper, which actually never
explicitly presents a procedure for matrix bandwidth recognition [CS05]. However, the paper
does define a bandwidth minimization algorithm that repeatedly calls a recognition
subroutine—this is precisely the logic we implement here. (We do, however, also implement
said minimization algorithm in
[`MatrixBandwidth.Minimization.Exact.CapraraSalazarGonzalez`](@ref).)

# References
- [CS05](@cite): A. Caprara and J.-J. Salazar-González. *Laying Out Sparse Graphs with
    Provably Minimum Bandwidth*. INFORMS Journal on Computing **17**, 356–73 (2005).
    https://doi.org/10.1287/ijoc.1040.0083.
"""
struct CapraraSalazarGonzalez <: AbstractDecider end

# push!(ALGORITHMS[:Recognition], CapraraSalazarGonzalez)

Base.summary(::CapraraSalazarGonzalez) = "Caprara–Salazar-González"

_requires_structural_symmetry(::CapraraSalazarGonzalez) = true

function _bool_bandwidth_k_ordering(
    A::AbstractMatrix{Bool}, k::Integer, ::CapraraSalazarGonzalez
)
    error("TODO: Not yet implemented")
    return nothing
end
