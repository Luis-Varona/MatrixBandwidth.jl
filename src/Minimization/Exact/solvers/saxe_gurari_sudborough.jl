# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SaxeGurariSudborough <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

[TODO: Write here]

# Supertype Hierarchy
`SaxeGurariSudborough` <: [`ExactSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

# Performance
Given an ``n×n`` input matrix ``A``, the Saxe–Gurari–Sudborough algorithm runs in
``O(nⁿ⁻¹)`` time:
- For each underlying "bandwidth ≤ k" check, we call the Saxe–Gurari–Sudborough recognition
  algorithm, which runs in ``O(nᵏ)`` time [GS84, p. 531]. (This is an improvement upon the
  original ``O(nᵏ⁺¹)`` Saxe algorithm [Sax80, p. 363].)
- The difference between the maximum possible bandwidth (``n - 1``) and our initial lower
    bound grows linearly in ``n``, so we run the underlying ``O(nᵏ)`` recognition algorithm
    ``O(n)`` times.
- Therefore, the overall time complexity is ``∑ₖ₌₀ⁿ⁻¹ nᵏ = O(nⁿ⁻¹)``.

# Examples
[TODO: Write here]

# Notes
Whereas most exact bandwidth minimization algorithms are technically factorial-time (with
respect to ``n``) in the worst case but practically always approximate exponential time
complexity in real life, the ``O(nⁿ⁻¹)`` upper bound on Saxe–Gurari–Sudborough is typically
a good representation of actual performance in most cases. Indeed, these other types of
algorithms tend to outperform Saxe–Gurari–Sudborough for larger ``n``, given that their
aggressive pruning strategies keep their effective search space very small in practice and
``O(2ⁿ)`` ⊂ ``O(nⁿ⁻¹)``.

# References
- [GS84](@cite): E. M. Gurari and I. H. Sudborough. *Improved dynamic programming algorithms
    for bandwidth minimization and the MinCut Linear Arrangement problem*. Journal of
    Algorithms **5**, 531–46 (1984). https://doi.org/10.1016/0196-6774(84)90006-3.
- [Sax80](@cite): J. B. Saxe. *Dynamic-Programming Algorithms for Recognizing
    Small-Bandwidth Graphs in Polynomial Time*. SIAM Journal on Algebraic and Discrete
    Methods **1**, 363–69 (1980). https://doi.org/10.1137/0601042.
"""
struct SaxeGurariSudborough <: ExactSolver end

push!(ALGORITHMS[:Minimization][:Exact], SaxeGurariSudborough)

Base.summary(::SaxeGurariSudborough) = "Saxe–Gurari–Sudborough"

_requires_structural_symmetry(::SaxeGurariSudborough) = true

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::SaxeGurariSudborough)
    components = _connected_components(A)
    #= Heuristically, it is likelier that smaller components have lower bandwidth, so we
    process them first to keep `k` low for as long as possible (since the complexity of each
    recognition subroutine is exponential in `k`). =#
    sort!(components; by=length)

    ordering = Vector{Int}(undef, size(A, 1))
    k = 1
    num_placed = 0

    for component in components
        submatrix = view(A, component, component)
        k = max(k, bandwidth_lower_bound(submatrix))
        component_ordering = Recognition._sgs_connected_ordering(submatrix, k)

        while isnothing(component_ordering)
            component_ordering = Recognition._sgs_connected_ordering(submatrix, k += 1)
        end

        ordering[(num_placed + 1):(num_placed += length(component))] .= component[component_ordering]
    end

    return ordering
end
