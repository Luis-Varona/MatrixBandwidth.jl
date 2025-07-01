# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    CapraraSalazarGonzalez <: ExactSolver <: AbstractSolver

TODO: Write here
"""
struct CapraraSalazarGonzalez <: ExactSolver end

Base.summary(::CapraraSalazarGonzalez) = "Caprara–Salazar-González algorithm"

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::CapraraSalazarGonzalez)
    # TODO: Implement
    return Int[]
end
