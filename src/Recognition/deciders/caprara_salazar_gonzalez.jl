# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CapraraSalazarGonzalez <: AbstractDecider

TODO: Write here
"""
struct CapraraSalazarGonzalez <: AbstractDecider end

Base.summary(::CapraraSalazarGonzalez) = "Caprara–Salazar-González"

_requires_symmetry(::CapraraSalazarGonzalez) = true

function _bool_bandwidth_k_ordering(
    A::AbstractMatrix{Bool}, k::Int, ::CapraraSalazarGonzalez
)
    # TODO: Implement
end
