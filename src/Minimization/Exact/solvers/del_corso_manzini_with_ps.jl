# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    DelCorsoManziniWithPS <: ExactSolver <: AbstractSolver

TODO: Write here
"""
struct DelCorsoManziniWithPS <: ExactSolver
    depth::Int

    DelCorsoManziniWithPS() = new(0)

    function DelCorsoManziniWithPS(depth::Int)
        if depth < 0
            throw(ArgumentError("MB-PS depth parameter must be positive, got $depth"))
        end

        return new(depth)
    end
end

Base.summary(::DelCorsoManziniWithPS) = "Del Corso–Manzini with perimeter search"

_requires_symmetry(::DelCorsoManziniWithPS) = true

function _optimal_mbps_depth(A::AbstractMatrix{Bool})
    n = size(A, 1)

    # I'm really not sure about this... let's check the original paper again
    if n <= 2
        return 1
    elseif n <= 10
        return 3
    elseif n <= 20
        return 4
    elseif n <= 50
        return 5
    elseif n <= 60
        return 7
    elseif n <= 80
        return 6
    else
        # TODO: Maybe implement sparsity-based heuristics?
    end
end

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, solver::DelCorsoManziniWithPS)
    depth = solver.depth

    if depth == 0
        depth = _optimal_mbps_depth(A)
    end

    # TODO: Implement

    return Int[] # Just a placeholder
end
