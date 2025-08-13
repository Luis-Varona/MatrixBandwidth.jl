# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    AbstractAlgorithm

Abstract base type for all matrix bandwidth minimization and recognition algorithms.

# Interface
Concrete subtypes of `AbstractAlgorithm` must implement the following methods:
- `Base.summary(::T) where {T<:AbstractAlgorithm}`: returns a `String` indicating the name
    of the algorithm (e.g., `"Gibbs–Poole–Stockmeyer"`).
- `_requires_symmetry(::T) where {T<:AbstractAlgorithm}`: returns a `Bool` indicating
    whether the algorithm requires the input matrix to be structurally symmetric.

Direct subtypes of `AbstractAlgorithm` must implement the following method:
- `_problem(::T) where {T<:AbstractAlgorithm}`: returns a `Symbol` indicating the
    matrix bandwidth problem tackled by the algorithm (e.g., `:minimization`).
"""
abstract type AbstractAlgorithm end

function Base.summary(::T) where {T<:AbstractAlgorithm}
    throw(NotImplementedError(summary, T, AbstractAlgorithm))
end

# Indicate whether the algorithm requires the input matrix to be structurally symmetric
function _requires_symmetry(::T) where {T<:AbstractAlgorithm}
    throw(NotImplementedError(_requires_symmetry, T, AbstractAlgorithm))
end

#= Indicate the matrix bandwidth problem tackled (e.g., `:minimization`). Each direct
subtype of `AbstractAlgorithm` must implement its own `_problem` method. =#
function _problem(::T) where {T<:AbstractAlgorithm}
    subtype = _find_direct_subtype(AbstractAlgorithm, T)
    throw(NotImplementedError(_problem, subtype, AbstractAlgorithm))
end

"""
    AbstractResult

Abstract base type for all matrix bandwidth problem results.

# Interface
Concrete subtypes of `AbstractResult` *must* implement parametric types
- `A<:AbstractAlgorithm`;
- `M<:AbstractMatrix{<:Number}`; and
- `O<:Union{Nothing,Vector{Int}}`,

alongside the following fields:
- `algorithm::A`: the algorithm used to investigate the bandwidth.
- `matrix::M`: the matrix whose bandwidth is investigated.
- `ordering::O`: the corresponding ordering of the rows and columns, if a relevant one is
    found; otherwise, `nothing`.
"""
abstract type AbstractResult end

Base.summary(res::AbstractResult) = summary(res.algorithm)

"""
    NotImplementedError{Nothing}(f, subtype, abstracttype)
    NotImplementedError{Symbol}(f, arg, subtype, abstracttype)

An exception indicating that a function lacks dispatch to handle a specific argument type.

Semantically, this differs from `MethodError` in that it connotes a developer-side failure
to implement a method rather than erroneous user input. Throughout this package, it is often
used to warn when an existing function with multiple dispatch on some abstract type is
called on a newly created subtype for which no method has been defined.

# Fields
- `f::Function`: the function called.
- `arg::Symbol`: the name of the argument with the unsupported type, if the function has
    multiple arguments. If the function has only one argument, this field should be set to
    `nothing`.
- `subtype::Type`: the type of the argument. May be the actual concrete type or some
    intermediate supertype. (For instance, if the relevant input has concrete type `A` with
    hierarchy `A <: B <: C` and the `abstracttype` field is `C`, then both `A` and `B` are
    perfectly valid choices for `subtype`.)
- `abstracttype::Type`: the abstract type under which the argument is meant to fall.

# Constructors
- `NotImplementedError(::Function, ::Type, ::Type)`: constructs a new `NotImplementedError`
    instance for a single-argument function. Throws an error if the second type is not
    abstract or the first type is not a subtype of the second.
- `NotImplementedError(::Function, ::Symbol, ::Type, ::Type)`: constructs a new
    `NotImplementedError` instance for a multi-argument function. Throws an error if the
    second type is not abstract or the first type is not a subtype of the second.

# Supertype Hierarchy
`NotImplementedError` <: `Exception`
"""
struct NotImplementedError{T<:Union{Nothing,Symbol}} <: Exception
    f::Function
    arg::T
    subtype::Type
    abstracttype::Type

    function NotImplementedError(f::Function, subtype::Type, abstracttype::Type)
        return NotImplementedError(f, nothing, subtype, abstracttype)
    end

    function NotImplementedError(
        f::Function, arg::T, subtype::Type, abstracttype::Type
    ) where {T<:Union{Nothing,Symbol}}
        if !isabstracttype(abstracttype)
            throw(ArgumentError("Expected an abstract type, got $abstracttype"))
        end

        if !(subtype <: abstracttype)
            throw(ArgumentError("Expected a subtype of $abstracttype, got $subtype"))
        end

        return new{T}(f, arg, subtype, abstracttype)
    end
end

function Base.showerror(io::IO, e::NotImplementedError{Nothing})
    print(
        io,
        """NotImplementedError with $(e.subtype):
        $(e.f) is not yet implemented for this subtype of $(e.abstracttype).
        Try defining method dispatch manually if this is a newly created subtype.""",
    )
    return nothing
end

function Base.showerror(io::IO, e::NotImplementedError{Symbol})
    print(
        io,
        """NotImplementedError with argument $(e.arg)::$(e.subtype):
        $(e.f) is not yet implemented for this subtype of $(e.abstracttype).
        Try defining method dispatch manually if this is a newly created subtype.""",
    )
    return nothing
end

"""
    RectangularMatrixError(A)

An exception indicating that the matrix `A` is not square.

Matrix bandwidth is only defined for square matrices, so this exception is raised when a
bandwidth minimization or recognition algorithm is called with a non-square input.

# Fields
- `A::AbstractMatrix{<:Number}`: the input matrix.
- `m::Int`: the number of rows of `A`.
- `n::Int`: the number of columns of `A`.

# Constructors
- `RectangularMatrixError(A::AbstractMatrix{<:Number})`: constructs a new
    `RectangularMatrixError` instance, automatically inferring `m` and `n` by calling
    `size(A)`.

# Supertype Hierarchy
`RectangularMatrixError` <: `Exception`
"""
struct RectangularMatrixError <: Exception
    A::AbstractMatrix{<:Number}
    m::Int
    n::Int

    function RectangularMatrixError(A::AbstractMatrix{<:Number})
        return new(A, size(A)...)
    end
end

function Base.showerror(io::IO, e::RectangularMatrixError)
    compact = get(io, :compact, true)::Bool
    limit = get(io, :limit, true)::Bool

    print(
        IOContext(io, :compact => compact, :limit => limit),
        "RectangularMatrixError with ",
        e.A,
    )
    print(
        io,
        """:
        Bandwidth is only defined for square matrices, got dimensions $(e.m)×$(e.n).""",
    )

    return nothing
end

"""
    StructuralAsymmetryError(A, algorithm)

An exception indicating that the matrix `A` is not structurally symmetric.

An `n×n` matrix ``A`` is *structurally symmetric* if ``A[i, j]`` is nonzero if and only if
``A[j, i]`` is nonzero for ``1 ≤ i, j ≤ n``. Many (albeit not all) matrix bandwidth
minimization and recognition algorithms assume structural symmetry, so this exception is
raised when one of these algorithms is called with a structurally asymmetric input.

# Fields
- `A::AbstractMatrix{<:Number}`: the input matrix.
- `algorithm::AbstractAlgorithm`: the algorithm that was called.
- `problem::Symbol`: the matrix bandwidth problem tackled by the algorithm (e.g.,
    `:minimization`).

# Constructors
- `StructuralAsymmetryError(A::AbstractMatrix{<:Number}, algorithm::AbstractAlgorithm)`:
    constructs a new `StructuralAsymmetryError` instance, automatically inferring `problem`
    by calling `_problem(algorithm)`.

# Supertype Hierarchy
`StructuralAsymmetryError` <: `Exception`

# Notes
As noted in the error message for `StructuralAsymmetryError` instances, users may want to
consider symmetrization techniques from [RS06] to minimize the bandwidth of structurally
asymmetric matrices. (A prominent one is to simply replace ``A[i, j]`` with ``1`` whenever
``A[i, j] = 0`` but ``A[j, i] ≠ 0``.) Of course, the reliability of minimization algorithms
is diminished after such a transformation, so users should proceed with caution nonetheless.

# References
- [RS06](@cite): J. K. Reid and J. A. Scott. *Reducing the Total Bandwidth of a Sparse
    Unsymmetric Matrix*. SIAM Journal on Matrix Analysis and Applications **28**, 805–21
    (2006). https://doi.org/10.1137/050629938.
"""
struct StructuralAsymmetryError <: Exception
    A::AbstractMatrix{<:Number}
    algorithm::AbstractAlgorithm
    problem::Symbol

    function StructuralAsymmetryError(
        A::AbstractMatrix{<:Number}, algorithm::AbstractAlgorithm
    )
        return new(A, algorithm, _problem(algorithm))
    end
end

function Base.showerror(io::IO, e::StructuralAsymmetryError)
    compact = get(io, :compact, true)::Bool
    limit = get(io, :limit, true)::Bool

    print(
        IOContext(io, :compact => compact, :limit => limit),
        "StructuralAsymmetryError with ",
        e.A,
    )
    print(
        io,
        """:
        $(Base.summary(e.algorithm)) can only handle bandwidth $(e.problem) for structurally symmetric matrices.
        You may want to consider symmetrization techniques from Reid and Scott (2006): https://doi.org/10.1137/050629938.""",
    )

    return nothing
end
