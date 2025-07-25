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
``k`` until a bandwidth-``k`` ordering is found [CSG05](@cite), with ``k`` initialized to
some lower bound on the minimum bandwidth of ``A`` up to symmetric permutation.

As noted above, the Caprara–Salazar-González algorithm requires structurally symmetric input
(that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for
``1 ≤ i, j ≤ n``).

# Performance
[TODO: Write here]

# Examples
[TODO: Write here]

# Notes
This algorithm is not the main one described in the original paper [CSG05](@cite), which
actually never explicitly presents a procedure for matrix bandwidth recognition. However,
the paper does define a bandwidth minimization algorithm that repeatedly calls a recognition
subroutine—this is precisely the logic we implement here. (We do, however, also implement
said minimization algorithm in
[`MatrixBandwidth.Minimization.Exact.CapraraSalazarGonzalez`](@ref).)
"""
struct CapraraSalazarGonzalez <: AbstractDecider end

Base.summary(::CapraraSalazarGonzalez) = "Caprara–Salazar-González"

_requires_symmetry(::CapraraSalazarGonzalez) = true

function _bool_bandwidth_k_ordering(
    A::AbstractMatrix{Bool}, k::Int, ::CapraraSalazarGonzalez
)
    error("TODO: Not yet implemented")
    return nothing
end
