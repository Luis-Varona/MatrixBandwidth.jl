# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SimulatedAnnealing <: MetaheuristicSolver <: AbstractSolver <: AbstractAlgorithm

[TODO: Write here]

# Supertype Hierarchy
`SimulatedAnnealing` <: [`MetaheuristicSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

[TODO: Write here]
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

_requires_symmetry(::SimulatedAnnealing) = false

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::SimulatedAnnealing)
    error("TODO: Not yet implemented")
    return nothing
end
