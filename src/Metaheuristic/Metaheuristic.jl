# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Metaheuristic

Metaheuristic solvers for matrix bandwidth minimization.

This submodule is part of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Metaheuristic

#! format: off
import ..AbstractSolver, ..NotImplementedError
import .._approach, .._bool_minimal_band_ordering
#! format: on

export SimulatedAnnealing, GeneticAlgorithm, GRASP

include("types.jl")

include("simulated_annealing.jl")
include("genetic_algorithm.jl")
include("grasp.jl")

end
