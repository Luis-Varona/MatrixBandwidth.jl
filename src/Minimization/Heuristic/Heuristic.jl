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
heuristic algorithms (such as reverse Cuthill–McKee) run in polynomial time.

Heuristic algorithms differ from metaheuristic ones in that they do not employ higher-level
iterative search frameworks (e.g., stochastic techniques) to survey the global search space
and escape local minima; instead, they rely on straightforward deterministic procedures to
find good solutions in a single pass.

The following heuristic algorithms are currently supported:
- Gibbs–Poole–Stockmeyer algorithm ([`GibbsPooleStockmeyer`](@ref))
- Cuthill–McKee algorithm ([`CuthillMcKee`](@ref))
- Reverse Cuthill–McKee algorithm ([`ReverseCuthillMcKee`](@ref))

This submodule is part of the `MatrixBandwidth.Minimization` submodule of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Heuristic

#! format: off
import ..AbstractSolver, ..NotImplementedError
import .._approach, .._assert_matrix_is_square, .._bool_minimal_band_ordering, .._symmetrize
#! format: on

using DataStructures: Queue, enqueue!, dequeue!

export CuthillMcKee, ReverseCuthillMcKee, GibbsPooleStockmeyer

include("utils.jl")
include("types.jl")

include("solvers/cuthill_mckee.jl")
include("solvers/reverse_cuthill_mckee.jl")
include("solvers/gibbs_poole_stockmeyer.jl")

const DEFAULT_SELECTOR = pseudo_peripheral_node

end
