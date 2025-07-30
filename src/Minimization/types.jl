# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    AbstractSolver <: AbstractAlgorithm

Abstract base type for all matrix bandwidth minimization solvers.

# Interface
As per the interface of supertype [`AbstractAlgorithm`](@ref), concrete subtypes of
`AbstractSolver` must implement the following methods:
- `Base.summary(::T) where {T<:AbstractSolver}`: returns a `String` indicating the name
    of the solver (e.g., `"Gibbs–Poole–Stockmeyer"`).
- `_requires_symmetry(::T) where {T<:AbstractSolver}`: returns a `Bool` indicating
    whether the solver requires the input matrix to be structurally symmetric.

Direct subtypes of `AbstractSolver` must implement the following method:
- `_approach(::T) where {T<:AbstractSolver}`: returns a `Symbol` indicating the
    category of solver (e.g., `:heuristic`).

# Supertype Hierarchy
`AbstractSolver` <: [`AbstractAlgorithm`](@ref)
"""
abstract type AbstractSolver <: AbstractAlgorithm end

_problem(::AbstractSolver) = :minimization

#= Indicate the category of solver (e.g., heuristic). Each direct subtype of
`AbstractSolver` must implement its own `_approach` method. =#
function _approach(::T) where {T<:AbstractSolver}
    subtype = _find_direct_subtype(AbstractSolver, T)
    throw(NotImplementedError(_approach, subtype, AbstractSolver))
end

"""
    MinimizationResult{A,M,O} <: AbstractResult

Output struct for matrix bandwidth minimization results.

# Fields
- `algorithm::A<:AbstractSolver`: the solver used to minimize the bandwidth.
- `matrix::M<:AbstractMatrix{<:Number}`: the original matrix whose bandwidth is minimized.
- `ordering::O<:Vector{Int}`: the (near-)optimal ordering of the rows and columns.
- `bandwidth::Int`: the minimized bandwidth of the matrix.
- `approach::Symbol`: the approach used by the solver. (Should be one of `:exact`,
    `:heuristic`, and `:metaheuristic`.)

# Constructors
- `MinimizationResult(algorithm, matrix, ordering, bandwidth)`: constructs a new
    `MinimizationResult` instance with the given fields. The `approach` field is
    automatically determined based on the algorithm type.

# Supertype Hierarchy
`MinimizationResult` <: [`AbstractResult`](@ref)
"""
struct MinimizationResult{A<:AbstractSolver,M<:AbstractMatrix{<:Number},O<:Vector{Int}} <:
       AbstractResult
    algorithm::A
    matrix::M
    ordering::O
    bandwidth::Int
    approach::Symbol # :exact, :heuristic, or :metaheuristic

    function MinimizationResult(
        algorithm::A, matrix::M, ordering::Vector{Int}, bandwidth::Int
    ) where {A<:AbstractSolver,M<:AbstractMatrix{<:Number}}
        return new{A,M,Vector{Int}}(
            algorithm, matrix, ordering, bandwidth, _approach(algorithm)
        )
    end
end

# The `Base.show` override here takes heavy inspiration from the `Optim.jl` package
function Base.show(io::IO, res::MinimizationResult)
    n = size(res.matrix, 1)

    println(io, "Results of Bandwidth Minimization Algorithm")
    println(io, " * Algorithm: $(summary(res.algorithm))")
    println(io, " * Approach: $(string(res.approach))")
    println(io, " * Minimum Bandwidth: $(res.bandwidth)")
    println(io, " * Original Bandwidth: $(bandwidth(res.matrix))")
    print(io, " * Matrix Size: $n×$n")

    return nothing
end
