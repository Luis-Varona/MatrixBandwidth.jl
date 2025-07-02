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
must implement its own `_approach` method for use in `BandMinResult` instantiation. =#
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
    BandMinResult

Output struct for matrix bandwidth minimization results.

# Fields
- `matrix::M`: the original matrix whose bandwidth is minimized.
- `bandwidth::Int`: the minimized bandwidth of the matrix.
- `ordering::Vector{Int}`: the (near-)optimal ordering of the rows and columns.
- `solver::S`: the algorithm used to minimize the bandwidth.
- `approach::Symbol`: the approach used by the solver. (Should be one of `:exact`,
    `:heuristic`, and `:metaheuristic`.)

# Constructors
- `BandMinResult(matrix, bandwidth, ordering, solver)`: constructs a new
    [`BandMinResult`](@ref) instance with the given fields. The `approach` field is
    automatically determined based on the solver type.
"""
struct BandMinResult{M<:AbstractMatrix{<:Number},S<:AbstractSolver}
    matrix::M
    bandwidth::Int
    ordering::Vector{Int}
    solver::S
    approach::Symbol # :exact, :heuristic, or :metaheuristic

    function BandMinResult(
        matrix::M, bandwidth::Int, ordering::Vector{Int}, solver::S
    ) where {M<:AbstractMatrix{<:Number},S<:AbstractSolver}
        return new{M,S}(matrix, bandwidth, ordering, solver, _approach(solver))
    end
end

Base.summary(res::BandMinResult) = summary(res.solver)

# The `Base.show` override here takes heavy inspiration from the `Optim.jl` package
function Base.show(io::IO, res::BandMinResult)
    n = size(res.matrix, 1)

    println(io, "Results of Bandwidth Minimization Algorithm")
    println(io, " * Algorithm: $(summary(res.solver))")
    println(io, " * Approach: $(string(res.approach))")
    println(io, " * Minimum Bandwidth: $(res.bandwidth)")
    println(io, " * Original Bandwidth: $(bandwidth(res.matrix))")
    print(io, " * Matrix Size: $n×$n")

    return nothing
end
