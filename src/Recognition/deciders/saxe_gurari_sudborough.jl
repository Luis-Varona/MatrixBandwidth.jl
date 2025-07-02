# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    SaxeGurariSudborough <: AbstractDecider

TODO: Write here
"""
struct SaxeGurariSudborough <: AbstractDecider end

Base.summary(::SaxeGurariSudborough) = "Saxe–Gurari–Sudborough algorithm"

function _bool_bandwidth_k_ordering(A::AbstractMatrix{Bool}, k::Int, ::SaxeGurariSudborough)
    # TODO: Implement
end
