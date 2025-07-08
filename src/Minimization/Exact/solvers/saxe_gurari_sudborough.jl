# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SaxeGurariSudborough <: ExactSolver <: AbstractSolver

TODO: Write here
"""
struct SaxeGurariSudborough <: ExactSolver end

Base.summary(::SaxeGurariSudborough) = "Saxe–Gurari–Sudborough"

_requires_symmetry(::SaxeGurariSudborough) = true

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::SaxeGurariSudborough)
    # TODO: Implement
    return Int[]
end
