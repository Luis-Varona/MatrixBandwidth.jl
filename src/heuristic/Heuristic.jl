# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Heuristic

Heuristic solvers for matrix bandwidth minimization.

This submodule is part of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Heuristic

#! format: off
import ..AbstractSolver, ..NotImplementedError, ..approach, .._sym_minimal_band_ordering
#! format: on

using DataStructures: Queue, enqueue!, dequeue!

export CuthillMcKee, ReverseCuthillMcKee

include("utils.jl")
include("types.jl")

include("cuthill_mckee.jl")
include("reverse_cuthill_mckee.jl")

const DEFAULT_SELECTOR = pseudo_peripheral_node

end
