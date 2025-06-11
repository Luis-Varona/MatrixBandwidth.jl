# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    NotImplementedError(f, arg, subtype, abstracttype)

An exception indicating that a function lacks dispatch to handle a specific argument type.

Semantically, this differs from `MethodError` in that it connotes a developer-side failure
to implement a method rather than erroneous user input. Throughout this package, it is often
used to warn when an existing function with multiple dispatch on some abstract type is
called on a newly created subtype for which no method has been defined.

# Fields
- `f::Function`: the function called.
- `arg::Symbol`: the name of the argument with the unsupported type.
- `subtype::Type`: the type of the argument. May be the actual concrete type or some
    intermediate supertype. (For instance, if the relevant input has concrete type `A` with
    hierarchy `A <: B <: C` and the `abstracttype` field is `C`, then both `A` and `B` are
    perfectly valid choices for `subtype`.)
- `abstracttype::Type`: the abstract type under which the argument is meant to fall.

# Constructors
- `NotImplementedError(::Function, ::Symbol, ::Type, ::Type)`: constructs a new
    `NotImplementedError` instance. Throws an error if the second type is not abstract or
    the first type is not a subtype of the second.
"""
struct NotImplementedError <: Exception
    f::Function
    arg::Symbol
    subtype::Type
    abstracttype::Type

    function NotImplementedError(
        f::Function, arg::Symbol, subtype::Type, abstracttype::Type
    )
        if !isabstracttype(abstracttype)
            throw(ArgumentError("Expected an abstract type, got $abstracttype"))
        end

        if !(subtype <: abstracttype)
            throw(ArgumentError("Expected a subtype of $abstracttype, got $subtype"))
        end

        return new(f, arg, subtype, abstracttype)
    end
end

function Base.showerror(io::IO, e::NotImplementedError)
    return print(
        io,
        """NotImplementedError with argument $(e.arg)::$(e.subtype):
        $(e.f) is not yet implemented for this subtype of $(e.abstracttype).
        Try defining method dispatch manually if this is a newly created subtype.""",
    )
end

"""
    random_sparse_banded_matrix(n, k; p=0.75, rng=default_rng()) -> Matrix{Float64}

Generate a random `n×n` matrix with bandwidth exactly `k` and sparse bands with density `p`.

All entries from this matrix will be from the interval `[0, 1]`. Entries up to the `k`-th
superdiagonal and down to the `k`-th subdiagonal are nonzero with probability `p`, and each
band has at least one nonzero entry to ensure that the bandwidth is precisely `k`.

# Arguments
- `n::Int`: the order of the matrix to generate. Must be positive.
- `k::Int`: the desired matrix bandwidth. Must satisfy `0 ≤ k < n`.

# Keyword Arguments
- `p::Float64=0.75`: the band density. Must satisfy `0 < p ≤ 1`.
- `rng::AbstractRNG=Random.default_rng()`: the random number generator to use. Defaults to
    `Random.default_rng()`.

# Returns
- `::Matrix{Float64}`: a random `n×n` matrix with bandwidth exactly `k` and sparse bands
    with density `p`.

# Examples
TODO: Write here

# Notes
Users of [`MatrixBandwidth`](@ref) may find this function useful when generating random test
data for whatever frameworks, algorithms, etc. they are implementing.
"""
function random_sparse_banded_matrix(
    n::Int, k::Int; p::Float64=0.75, rng::AbstractRNG=Random.default_rng()
)
    if n <= 0
        throw(ArgumentError("Matrix order must be positive, got $n"))
    end

    if k < 0
        throw(ArgumentError("Matrix bandwidth must be non-negative, got $k"))
    end

    if n <= k
        throw(ArgumentError("Matrix order must be greater than bandwidth, got $n and $k"))
    end

    if !(0 < p <= 1)
        throw(ArgumentError("Band density must be in (0, 1], got $p"))
    end

    A = zeros(Float64, n, n)

    #= Ensure that each band has at least one nonzero entry. (The `is_upper` flag indicates
    whether we are filling a superdiagonal or subdiagonal.) =#
    for offset in 0:k, is_upper in (true, false)
        rows = ((1 + offset):n) .- is_upper * offset
        cols = (1:(n - offset)) .+ is_upper * offset
        idx = rand(rng, 1:length(rows))
        A[rows[idx], cols[idx]] = rand(rng)
    end

    for i in 1:n, j in max(1, i - k):(min(n, i + k))
        if rand(rng) < p
            A[i, j] = rand(rng)
        end
    end

    return A
end

# Validate that some matrix is square before attempting to minimize/compute its bandwidth
function _assert_matrix_is_square(A::AbstractMatrix{T}) where {T<:Number}
    if !allequal(size(A))
        throw(ArgumentError("Matrix bandwidth is not defined for non-square matrices"))
    end
end
