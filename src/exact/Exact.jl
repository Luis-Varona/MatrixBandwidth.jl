# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Exact

Exact solvers for matrix bandwidth minimization.

This submodule is part of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Exact

#! format: off
import ..AbstractSolver, ..NotImplementedError, .._approach, .._bool_minimal_band_ordering
#! format: on

export MBID, MBPS

include("types.jl")

include("mbid.jl")
include("mbps.jl")

end
