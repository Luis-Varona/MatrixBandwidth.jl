# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    GibbsPooleStockmeyer <: HeuristicSolver <: AbstractSolver <: AbstractAlgorithm

[TODO: Write here]

As noted above, the Gibbs–Poole–Stockmeyer algorithm requires structurally symmetric input
(that is, ``A[i, j]`` must be nonzero if and only if ``A[j, i]`` is nonzero for
``1 ≤ i, j ≤ n``).

# Fields
- `node_selector::Function`: a function that selects a node from some connected component of
    the input matrix from which to start the breadth-first search. If no custom heuristic is
    specified, this field defaults to [`pseudo_peripheral_node`](@ref), which picks a node
    "farthest" from the others in the component (not necessarily the lowest-degree node).

# Performance
Given an ``n×n`` input matrix ``A``, the Cuthill–McKee algorithm runs in ``O(n²)`` time.

[TODO: Write here]

# Examples
[TODO: Write here]

# Notes
Note that the `node_selector` field must be of the form
`(A::AbstractMatrix{Bool}) -> Integer` (i.e., it must take in an boolean matrix and return
an integer). If this is not the case, an `ArgumentError` is thrown upon construction.
"""
struct GibbsPooleStockmeyer <: HeuristicSolver
    node_selector::Function

    function GibbsPooleStockmeyer(node_selector::Function=DEFAULT_SELECTOR)
        _assert_valid_node_selector(node_selector)
        return new(node_selector)
    end
end

Base.summary(::GibbsPooleStockmeyer) = "Gibbs–Poole–Stockmeyer"

_requires_symmetry(::GibbsPooleStockmeyer) = true

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::GibbsPooleStockmeyer)
    node_selector = solver.node_selector
    components = _connected_components(A)

    component_orderings = Iterators.map(
        component -> _gps_connected_ordering(view(A, component, component), node_selector),
        components,
    )

    return collect(
        Iterators.flatmap(
            ((component, component_ordering),) -> component[component_ordering],
            zip(components, component_orderings),
        ),
    )
end

# Gibbs–Poole–Stockmeyer searches each connected component independently.
function _gps_connected_ordering(A::AbstractMatrix{Bool}, node_selector::Function)
    error("TODO: Not yet implemented")
    return nothing
end
