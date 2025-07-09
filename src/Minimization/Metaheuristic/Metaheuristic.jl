# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Minimization.Metaheuristic

Metaheuristic solvers for matrix bandwidth minimization.

Metaheuristic methods are those which [TODO: Write here]

The following metaheuristic algorithms are currently supported:
- Greedy randomized adaptive search procedure (GRASP) ([`GRASP`](@ref))
- Simulated annealing ([`SimulatedAnnealing`](@ref))
- Genetic algorithm ([`GeneticAlgorithm`](@ref))

This submodule is part of the `MatrixBandwidth.Minimization` submodule of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Metaheuristic

#! format: off
import ..AbstractSolver
import ..NotImplementedError
import .._requires_symmetry
import .._approach, .._bool_minimal_band_ordering
#! format: on

export SimulatedAnnealing, GeneticAlgorithm, GRASP

include("types.jl")

include("solvers/simulated_annealing.jl")
include("solvers/genetic_algorithm.jl")
include("solvers/grasp.jl")

end
