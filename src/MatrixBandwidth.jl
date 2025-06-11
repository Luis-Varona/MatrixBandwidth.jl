# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth

Exact, heuristic, and metaheuristic algorithms for matrix bandwidth minimization in Julia.

The following algorithms are currently supported:
- **Exact**
    - [`MBID`](@ref): Minimum bandwidth by iterative deepening (MB-ID)
    - [`MBPS`](@ref): Minimum bandwidth by perimeter search (MB-PS)
- **Heuristic**
    - [`CuthillMcKee`](@ref): Cuthill–McKee algorithm
    - [`ReverseCuthillMcKee`](@ref): Reverse Cuthill–McKee algorithm
- **Metaheuristic**
    - [`SimulatedAnnealing`](@ref): Simulated annealing
    - [`GeneticAlgorithm`](@ref): Genetic algorithm
    - [`GRASP`](@ref): Greedy randomized adaptive search procedure (GRASP)

[Full documentation](https://Luis-Varona.github.io/MatrixBandwidth.jl/dev/) is available for
the latest development version of this package.
"""
module MatrixBandwidth

using Random

include("utils.jl")
include("types.jl")
include("minimize_bandwidth.jl")

include("exact/Exact.jl")
include("heuristic/Heuristic.jl")
include("metaheuristic/Metaheuristic.jl")

using .Exact, .Heuristic, .Metaheuristic

# The output struct, main minimization function, and raw bandwidth function
export BandwidthResult, minimize_bandwidth, bandwidth_unpermuted
export MBID, MBPS # Exact solvers
export CuthillMcKee, ReverseCuthillMcKee # Heuristic solvers
export SimulatedAnnealing, GeneticAlgorithm, GRASP # Metaheuristic solvers

const DEFAULT_SOLVER = ReverseCuthillMcKee()

end
