# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    AbstractSolver

Abstract base type for all matrix bandwidth minimization solvers.
"""
abstract type AbstractSolver end

#= Indicate the category of solver (e.g., heuristic). Each concrete `AbstractSolver` subtype
must implement its own `_approach` method for use in `BandwidthResult` instantiation. =#
function _approach(::T) where {T<:AbstractSolver}
    S = supertype(T)

    if S === AbstractSolver
        subtype = T
    else
        subtype = S
    end

    throw(NotImplementedError(_approach, :solver, subtype, AbstractSolver))
end

"""
    BandwidthResult

Output struct for matrix bandwidth minimization results.

# Fields
- `matrix::M`: the original matrix whose bandwidth is minimized.
- `bandwidth::Int`: the minimized bandwidth of the matrix.
- `ordering::Vector{Int}`: the (near-)optimal ordering of the rows and columns.
- `solver::S`: the algorithm used to minimize the bandwidth.
- `approach::Symbol`: the approach used by the solver. (Should be one of `:exact`,
    `:heuristic`, and `:metaheuristic`.)

# Constructors
- `BandwidthResult(matrix, bandwidth, ordering, solver)`: constructs a new `BandwidthResult`
    instance with the given fields. The `approach` field is automatically determined based
    on the solver type.
"""
struct BandwidthResult{M<:AbstractMatrix{<:Number},S<:AbstractSolver}
    matrix::M
    bandwidth::Int
    ordering::Vector{Int}
    solver::S
    approach::Symbol # :exact, :heuristic, or :metaheuristic

    function BandwidthResult(
        matrix::M, bandwidth::Int, ordering::Vector{Int}, solver::S
    ) where {M<:AbstractMatrix{<:Number},S<:AbstractSolver}
        return new{M,S}(matrix, bandwidth, ordering, solver, _approach(solver))
    end
end

Base.summary(res::BandwidthResult) = summary(res.solver)

# TODO: Override `Base.show(io::IO, res::BandwidthResult)`
