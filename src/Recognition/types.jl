# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    AbstractDecider

Abstract base type for all matrix bandwidth recognition deciders.
"""
abstract type AbstractDecider end

"""
    BandRecogResult

Output struct for matrix bandwidth recognition results.

# Fields
- `matrix::M`: the original matrix whose bandwidth is to be tested.
- `k::Int`: the threshold bandwidth against which to test.
- `bandwidth_k_ordering::Union{Nothing,Vector{Int}}`: an ordering of the rows and columns of
    `matrix` inducing a bandwidth at most `k`, if one exists; otherwise, `nothing`.
- `decider::D`: the algorithm used to test the bandwidth.
"""
struct BandRecogResult{
    M<:AbstractMatrix{<:Number},B<:Union{Nothing,Vector{Int}},D<:AbstractDecider
}
    matrix::M
    k::Int
    bandwidth_k_ordering::Union{Nothing,Vector{Int}}
    decider::D
end

Base.summary(res::BandRecogResult) = summary(res.decider)

# The `Base.show` override here takes heavy inspiration from the `Optim.jl` package
function Base.show(io::IO, res::BandRecogResult)
    n = size(res.matrix, 1)

    if isnothing(res.bandwidth_k_ordering)
        has_ordering = "no"
    else
        has_ordering = "yes"
    end

    println(io, "Results of Bandwidth Recognition Algorithm")
    println(io, " * Algorithm: $(summary(res.decider))")
    println(io, " * k: $(res.k)")
    println(io, " * Has Bandwidth ≤ k Ordering: $has_ordering")
    println(io, " * Original Bandwidth: $(bandwidth(res.matrix))")
    print(io, " * Matrix Size: $n×$n")

    return nothing
end
