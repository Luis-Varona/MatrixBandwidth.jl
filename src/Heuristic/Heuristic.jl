# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Heuristic

Heuristic solvers for matrix bandwidth minimization.

Heuristic methods are those which aim to produce near-optimal solutions in a more performant
manner than exact methods. While precise bandwidth minimization is NP-complete, many
heuristic algorithms (such as reverse Cuthill–McKee) run in polynomial time.

Heuristic algorithms differ from metaheuristic ones in that they do not employ higher-level
iterative search frameworks (e.g., stochastic techniques) to survey the global search space
and escape local minima; instead, they rely on straightforward deterministic procedures to
find good solutions in a single pass.

The following heuristic algorithms are currently supported:
- [`CuthillMcKee`](@ref): Cuthill–McKee algorithm
- [`ReverseCuthillMcKee`](@ref): Reverse Cuthill–McKee algorithm

This submodule is part of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Heuristic

#! format: off
import ..AbstractSolver, ..NotImplementedError
import .._approach, .._bool_minimal_band_ordering, .._symmetrize
#! format: on

using DataStructures: Queue, enqueue!, dequeue!

export CuthillMcKee, ReverseCuthillMcKee

include("utils.jl")
include("types.jl")

include("cuthill_mckee.jl")
include("reverse_cuthill_mckee.jl")

const DEFAULT_SELECTOR = pseudo_peripheral_node

end
