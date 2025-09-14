# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Minimization.Metaheuristic

Metaheuristic solvers for matrix bandwidth minimization.

Metaheuristic methods are those which employ higher-level iterative search frameworks such
as stochastic techniques or nature-inspired processes to survey the global search space and
escape local minima. Unlike heuristic methods—which follow fixed deterministic
procedures—metaheuristics adaptively refine candidate solutions over multiple iterations.
Although metaheuristic approaches are often slower than heuristic ones (but certainly still
faster than exact ones), they shine in complex cases where the latter may get trapped in
poor-quality local minima.

The following metaheuristic algorithms are currently supported:
- Greedy randomized adaptive search procedure (GRASP) ([`GRASP`](@ref))
- Simulated annealing ([`SimulatedAnnealing`](@ref))
- Genetic algorithm ([`GeneticAlgorithm`](@ref))

This submodule is part of the `MatrixBandwidth.Minimization` submodule of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Metaheuristic

using MatrixBandwidth
using MatrixBandwidth: NotImplementedError
using MatrixBandwidth: _requires_structural_symmetry

using MatrixBandwidth.Minimization
using MatrixBandwidth.Minimization: _approach, _bool_minimal_band_ordering

export
    # Types
    MetaheuristicSolver,

    # Metaheuristic solvers
    GRASP,
    SimulatedAnnealing,
    GeneticAlgorithm

MatrixBandwidth.ALGORITHMS[:Minimization][:Metaheuristic] = []

include("types.jl")

include("solvers/simulated_annealing.jl")
include("solvers/genetic_algorithm.jl")
include("solvers/grasp.jl")

end
