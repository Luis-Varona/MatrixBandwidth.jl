# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SimulatedAnnealing <: MetaheuristicSolver <: AbstractSolver

TODO: Write here
"""
struct SimulatedAnnealing <: MetaheuristicSolver
    initial_temp::Float64
    cooling_rate::Float64
    max_iterations::Int
    max_no_improve::Int
    seed::Union{Nothing,Int}

    # TODO: Make constructor with default values
end

Base.summary(::SimulatedAnnealing) = "Simulated annealing"

function _minimize_bandwidth_safe(A::AbstractMatrix{<:Bool}, solver::SimulatedAnnealing)
    # TODO: Implement
end
