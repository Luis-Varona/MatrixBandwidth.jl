# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    GRASP <: MetaheuristicSolver <: AbstractSolver <: AbstractAlgorithm

[TODO: Write here]

# Supertype Hierarchy
`GRASP` <: [`MetaheuristicSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

[TODO: Write here]
"""
struct GRASP <: MetaheuristicSolver
    # TODO: Define fields and constructor (for default values)
end

# push!(ALGORITHMS[:Minimization][:Metaheuristic], GRASP)

Base.summary(::GRASP) = "Greedy randomized adaptive search procedure (GRASP)"

_requires_symmetry(::GRASP) = false

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, Solver::GRASP)
    error("TODO: Not yet implemented")
    return nothing
end
