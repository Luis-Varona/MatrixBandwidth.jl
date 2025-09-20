# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    AntColony <: MetaheuristicSolver <: AbstractSolver <: AbstractAlgorithm

[TODO: Write here]

# Supertype Hierarchy
`AntColony` <: [`MetaheuristicSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

[TODO: Write here]
"""
struct AntColony <: MetaheuristicSolver
    # TODO: Define fields and constructor (for default values)
end

# push!(MatrixBandwidth.ALGORITHMS[:Minimization][:Metaheuristic], AntColony)

Base.summary(::AntColony) = "Ant colony optimization"

MatrixBandwidth._requires_structural_symmetry(::AntColony) = false

function Minimization._minimize_bandwidth_impl(A::AbstractMatrix{Bool}, solver::AntColony)
    error("TODO: Not yet implemented")
    return nothing
end
