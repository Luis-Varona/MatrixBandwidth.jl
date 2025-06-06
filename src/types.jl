# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    AbstractSolver

TODO: Write here
"""
abstract type AbstractSolver end

"""
    ExactSolver <: AbstractSolver

TODO: Write here
"""
abstract type ExactSolver <: AbstractSolver end

"""
    HeuristicSolver <: AbstractSolver

TODO: Write here
"""
abstract type HeuristicSolver <: AbstractSolver end

"""
    MetaheuristicSolver <: AbstractSolver

TODO: Write here
"""
abstract type MetaheuristicSolver <: AbstractSolver end

"""
    approach(solver::AbstractSolver) -> Symbol

TODO: Write here
"""
approach(::ExactSolver) = :exact
approach(::HeuristicSolver) = :heuristic
approach(::MetaheuristicSolver) = :metaheuristic

"""
    BandwidthResult

TODO: Write here
"""
struct BandwidthResult{M<:AbstractMatrix{<:Number},S<:AbstractSolver}
    matrix::M
    bandwidth::Int
    perm::Vector{Int}
    solver::S
    approach::Symbol # :exact, :heuristic, or :metaheuristic

    function BandwidthResult(
        matrix::M, bandwidth::Int, perm::Vector{Int}, solver::S
    ) where {M<:AbstractMatrix{<:Number},S<:AbstractSolver}
        return new{M,S}(matrix, bandwidth, perm, solver, approach(solver))
    end
end

Base.summary(res::BandwidthResult) = summary(res.solver)

# TODO: Override `Base.show(io::IO, res::BandwidthResult)`
