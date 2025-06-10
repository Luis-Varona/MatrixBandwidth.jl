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

include("utils.jl")
include("types.jl")
include("minimize_bandwidth.jl")

include("exact/Exact.jl")
include("heuristic/Heuristic.jl")
include("metaheuristic/Metaheuristic.jl")

using .Exact, .Heuristic, .Metaheuristic

export minimize_bandwidth # The main function
export MBID, MBPS # Exact solvers
export CuthillMcKee, ReverseCuthillMcKee # Heuristic solvers
export SimulatedAnnealing, GeneticAlgorithm, GRASP # Metaheuristic solvers

const DEFAULT_SOLVER = ReverseCuthillMcKee()

end
