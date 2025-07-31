# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SaxeGurariSudborough <: AbstractDecider <: AbstractAlgorithm

[TODO: Write here]

# Supertype Hierarchy
`SaxeGurariSudborough` <: [`AbstractDecider`](@ref) <: [`AbstractAlgorithm`](@ref)

[TODO: Write here]
"""
struct SaxeGurariSudborough <: AbstractDecider end

# push!(ALGORITHMS[:Recognition], SaxeGurariSudborough)

Base.summary(::SaxeGurariSudborough) = "Saxe–Gurari–Sudborough"

function _bool_bandwidth_k_ordering(
    A::AbstractMatrix{Bool}, k::Integer, ::SaxeGurariSudborough
)
    error("TODO: Not yet implemented")
    return nothing
end
