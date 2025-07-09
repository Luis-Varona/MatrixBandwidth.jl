# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    GibbsPooleStockmeyer <: HeuristicSolver <: AbstractSolver

TODO: Write here. Do we need a node selector?
"""
struct GibbsPooleStockmeyer <: HeuristicSolver end

Base.summary(::GibbsPooleStockmeyer) = "Gibbs–Poole–Stockmeyer"

_requires_symmetry(::GibbsPooleStockmeyer) = true

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::GibbsPooleStockmeyer)
    # TODO: Implement

    return Int[]
end
