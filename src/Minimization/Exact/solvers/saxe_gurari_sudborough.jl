# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SaxeGurariSudborough <: ExactSolver <: AbstractSolver <: AbstractAlgorithm

[TODO: Write here]

# Supertype Hierarchy
`SaxeGurariSudborough` <: [`ExactSolver`](@ref) <: [`AbstractSolver`](@ref) <: [`MatrixBandwidth.AbstractAlgorithm`](@ref)

[TODO: Write here]
"""
struct SaxeGurariSudborough <: ExactSolver end

# push!(ALGORITHMS[:Minimization][:Exact], SaxeGurariSudborough)

Base.summary(::SaxeGurariSudborough) = "Saxe–Gurari–Sudborough"

_requires_symmetry(::SaxeGurariSudborough) = true

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::SaxeGurariSudborough)
    components = _connected_components(A)
    sort!(components; by=length)

    ordering = Vector{Int}(undef, size(A, 1))
    k = bandwidth_lower_bound(A)
    num_placed = 0

    for component in components
        submatrix = view(A, component, component)
        component_ordering = Recognition._sgs_connected_ordering(submatrix, k)

        while isnothing(component_ordering)
            component_ordering = Recognition._sgs_connected_ordering(submatrix, k += 1)
        end

        ordering[(num_placed + 1):(num_placed += length(component))] .= component[component_ordering]
    end

    return ordering
end
