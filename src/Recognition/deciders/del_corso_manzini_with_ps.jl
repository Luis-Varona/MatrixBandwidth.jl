# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    DelCorsoManziniWithPS <: AbstractDecider

TODO: Write here
"""
struct DelCorsoManziniWithPS <: AbstractDecider
    depth::Int

    DelCorsoManziniWithPS() = new(0)

    function DelCorsoManziniWithPS(depth::Int)
        if depth < 0
            throw(ArgumentError("MB-PS depth parameter must be positive, got $depth"))
        end

        return new(depth)
    end
end

Base.summary(::DelCorsoManziniWithPS) = "Del Corso–Manzini algorithm with perimeter search"

# TODO: Implement the rest. Might need to extract stuff out from the minimization version.

function _bool_bandwidth_k_ordering(
    A::AbstractMatrix{Bool}, k::Int, decider::DelCorsoManziniWithPS
)
    # TODO: Implement
end
