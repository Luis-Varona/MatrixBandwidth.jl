# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TabuSearch <: MetaheuristicSolver <: AbstractSolver <: AbstractAlgorithm

[TODO: Write here]

# Supertype Hierarchy
`TabuSearch` <: [`MetaheuristicSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

[TODO: Write here]
"""
struct TabuSearch <: MetaheuristicSolver
    # TODO: Define fields and constructor (for default values)
end

# push!(MatrixBandwidth.ALGORITHMS[:Minimization][:Metaheuristic], TabuSearch)

Base.summary(::TabuSearch) = "Tabu search"

MatrixBandwidth._requires_structural_symmetry(::TabuSearch) = false

function Minimization._minimize_bandwidth_impl(A::AbstractMatrix{Bool}, solver::TabuSearch)
    error("TODO: Not yet implemented")
    return nothing
end
