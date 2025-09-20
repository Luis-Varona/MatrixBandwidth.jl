# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Minimization.Heuristic

Heuristic solvers for matrix bandwidth minimization.

Heuristic methods are those which aim to produce near-optimal solutions in a more performant
manner than exact methods. While precise bandwidth minimization is NP-complete, many
heuristic algorithms (such as Gibbs–Poole–Stockmeyer) run in polynomial time.

Heuristic algorithms differ from metaheuristic ones in that they do not employ higher-level
iterative search frameworks (e.g., stochastic techniques) to survey the global search space
and escape local minima; instead, they rely on straightforward deterministic procedures to
find good solutions in a single pass.

The following heuristic algorithms are currently supported:
- Gibbs–Poole–Stockmeyer ([`GibbsPooleStockmeyer`](@ref))
- Cuthill–McKee ([`CuthillMcKee`](@ref))
- Reverse Cuthill–McKee ([`ReverseCuthillMcKee`](@ref))

This submodule is part of the `MatrixBandwidth.Minimization` submodule of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Heuristic

using MatrixBandwidth
using MatrixBandwidth: NotImplementedError, StructuralAsymmetryError
using MatrixBandwidth: connected_components
using MatrixBandwidth: _requires_structural_symmetry

using MatrixBandwidth.Minimization
using MatrixBandwidth.Minimization: _approach, _minimize_bandwidth_impl

using DataStructures: Queue

export
    # Types
    HeuristicSolver,

    # Heuristic solvers
    GibbsPooleStockmeyer,
    CuthillMcKee,
    ReverseCuthillMcKee

MatrixBandwidth.ALGORITHMS[:Minimization][:Heuristic] = []

include("node_finders.jl")
include("types.jl")

include("solvers/gibbs_poole_stockmeyer.jl")
include("solvers/cuthill_mckee.jl") # Defines both `CuthillMcKee` and `ReverseCuthillMcKee`

const DEFAULT_NODE_FINDER = bi_criteria_node_finder

end
