# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    GibbsPooleStockmeyer <: HeuristicSolver <: AbstractSolver

TODO: Write here
"""
struct GibbsPooleStockmeyer <: HeuristicSolver end

Base.summary(::GibbsPooleStockmeyer) = "Gibbs–Poole–Stockmeyer algorithm"

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::GibbsPooleStockmeyer)
    # TODO: Implement
    return Int[]
end
