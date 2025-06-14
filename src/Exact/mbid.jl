# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MBID <: ExactSolver <: AbstractSolver

TODO: Write here
"""
struct MBID <: ExactSolver end

Base.summary(::MBID) = "Matrix bandwidth by iterative deepening"

function _bool_minimal_band_ordering(A::AbstractMatrix{Bool}, ::MBID)
    n = size(A, 1)
    max_degree = maximum(sum(A, dims=1))
    bandwidth = ceil(Int, max_degree / 2)
    found = false

    while !found
        unselected = Set(1:n)
        adjacency_list = Set{Int}()
    end
end

function _add_node!(unselected::Set{Int}, adjacency_list::Set{Int}, i::Int)
end
