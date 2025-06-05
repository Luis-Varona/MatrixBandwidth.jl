# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth

Exact, heuristic, and metaheuristic algorithms for matrix bandwidth minimization in Julia.

[Full documentation](https://Luis-Varona.github.io/MatrixBandwidth.jl/dev/) is available for
the latest development version of this package.
"""
module MatrixBandwidth

# TODO: Exports go here

include("minimize_bandwidth.jl")

include("exact/mbid.jl")
include("exact/mbps.jl")

include("heuristic/cuthill_mckee.jl")
include("heuristic/reverse_cuthill_mckee.jl")

include("metaheuristic/simulated_annealing.jl")
include("metaheuristic/genetic_algorithm.jl")
include("metaheuristic/grasp.jl")

end
