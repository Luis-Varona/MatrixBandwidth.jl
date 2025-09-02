# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    AbstractDecider <: AbstractAlgorithm

Abstract base type for all matrix bandwidth recognition deciders.

# Interface
As per the interface of supertype [`AbstractAlgorithm`](@ref), concrete subtypes of
`AbstractDecider` must implement the following methods:
- `Base.summary(::T) where {T<:AbstractDecider}`: returns a `String` indicating the name
    of the decider (e.g., `"Caprara–Salazar-González"`).
- `_requires_structural_symmetry(::T) where {T<:AbstractDecider}`: returns a `Bool`
    indicating whether the decider requires the input matrix to be structurally symmetric.

# Supertype Hierarchy
`AbstractDecider` <: [`AbstractAlgorithm`](@ref)
"""
abstract type AbstractDecider <: AbstractAlgorithm end

_problem(::AbstractDecider) = :recognition

"""
    RecognitionResult{A,M,O} <: AbstractResult

Output struct for matrix bandwidth recognition results.

# Fields
- `algorithm::A<:AbstractDecider`: the decider used to test the bandwidth.
- `matrix::M<:AbstractMatrix{<:Number}`: the original matrix whose bandwidth is tested.
- `ordering::O<:Union{Nothing,Vector{Int}}`: an ordering of the rows and columns of `matrix`
    inducing a bandwidth at most `k`, if such an ordering exists; otherwise, `nothing`.
- `k::Integer`: the threshold bandwidth against which to test.
- `has_ordering::Bool`: whether the matrix has an ordering inducing a bandwidth at most `k`.
    (This is `true` if and only if `ordering` is not `nothing`.)

# Constructors
- `RecognitionResult(decider, matrix, ordering, k)`: constructs a new `RecognitionResult`
    instance with the given fields. The `has_ordering` field is automatically determined
    based on whether `ordering` is `nothing` or a `Vector{Int}`.

# Supertype Hierarchy
`RecognitionResult` <: [`AbstractResult`](@ref)
"""
struct RecognitionResult{
    A<:AbstractDecider,M<:AbstractMatrix{<:Number},O<:Union{Nothing,Vector{Int}}
} <: AbstractResult
    algorithm::A
    matrix::M
    ordering::O
    k::Integer
    has_ordering::Bool

    function RecognitionResult(
        algorithm::A, matrix::M, ::Nothing, k::Integer
    ) where {A<:AbstractDecider,M<:AbstractMatrix{<:Number}}
        return new{A,M,Nothing}(algorithm, matrix, nothing, k, false)
    end

    function RecognitionResult(
        algorithm::A, matrix::M, ordering::Vector{Int}, k::Integer
    ) where {A<:AbstractDecider,M<:AbstractMatrix{<:Number}}
        return new{A,M,Vector{Int}}(algorithm, matrix, ordering, k, true)
    end
end

# The `Base.show` override here takes heavy inspiration from the Optim.jl package
function Base.show(io::IO, res::RecognitionResult)
    n = size(res.matrix, 1)

    println(io, "Results of Bandwidth Recognition Algorithm")
    println(io, " * Algorithm: $(summary(res.algorithm))")
    println(io, " * Bandwidth Threshold k: $(res.k)")
    println(io, " * Has Bandwidth ≤ k Ordering: $(res.has_ordering)")
    println(io, " * Original Bandwidth: $(bandwidth(res.matrix))")
    print(io, " * Matrix Size: $n×$n")

    return nothing
end
