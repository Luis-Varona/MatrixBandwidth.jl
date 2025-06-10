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
    approach(solver::AbstractSolver) -> Symbol

TODO: Write here
"""
function approach(::T) where {T<:AbstractSolver}
    S = supertype(T)

    if S === AbstractSolver
        subtype = T
    else
        subtype = S
    end

    throw(NotImplementedError(approach, :solver, subtype, AbstractSolver))
end

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
