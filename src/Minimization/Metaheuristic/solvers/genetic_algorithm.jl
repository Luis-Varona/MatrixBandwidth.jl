# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    GeneticAlgorithm <: MetaheuristicSolver <: AbstractSolver <: AbstractAlgorithm

[TODO: Write here]

# Supertype Hierarchy
`GeneticAlgorithm` <: [`MetaheuristicSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

[TODO: Write here]
"""
struct GeneticAlgorithm <: MetaheuristicSolver
    # TODO: Define fields and constructor (for default values)
end

# push!(MatrixBandwidth.ALGORITHMS[:Minimization][:Metaheuristic], GeneticAlgorithm)

Base.summary(::GeneticAlgorithm) = "Genetic algorithm"

MatrixBandwidth._requires_structural_symmetry(::GeneticAlgorithm) = false

function Minimization._bool_minimal_band_ordering(
    A::AbstractMatrix{Bool}, solver::GeneticAlgorithm
)
    error("TODO: Not yet implemented")
    return nothing
end
