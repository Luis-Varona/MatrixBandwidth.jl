# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Exact

Exact solvers for matrix bandwidth minimization.

Exact methods are those which guarantee an optimal ordering producing the true minimum
bandwidth of a matrix. Since bandwidth minimization is an NP-complete problem, existing
exact algorithms are, at best, exponential in time complexity—much worse than many
polynomial-time heuristic approaches (e.g., reverse Cuthill–McKee). Such methods, therefore,
are not feasible for large matrices, but they remain useful when precise solutions are
required for small-to-medium-sized inputs (say, up to ``100×100``).

The following exact algorithms are currently supported:
- [`MBID`](@ref): Minimum bandwidth by iterative deepening (MB-ID)
- [`MBPS`](@ref): Minimum bandwidth by perimeter search (MB-PS)

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
