# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth

Fast algorithms for matrix bandwidth minimization and recognition, written in Julia.

Reordering the rows and columns of a matrix to reduce its bandwidth has many practical
applications in engineering and scientific computing. It is a common preprocessing step used
to improve performance when solving linear systems, approximating partial differential
equations, optimizing circuit layout, and more [Maf14, p. 184].

Recall that the *bandwidth* of an ``n×n`` matrix ``A`` is the minimum non-negative integer
``k ∈ \\{0, 1, …, n - 1\\}`` such that ``A[i, j] = 0`` whenever ``|i - j| > k``.
Equivalently, ``A`` has bandwidth *at most* ``k`` if all entries above the ``k``ᵗʰ
superdiagonal and below the ``k``ᵗʰ subdiagonal are zero, and ``A`` has bandwidth *at least*
``k`` if there exists any nonzero entry in the ``k``ᵗʰ superdiagonal or subdiagonal.

The *matrix bandwidth minimization problem* involves finding a permutation matrix ``P`` such
that the bandwidth of ``PAPᵀ`` is minimized; this is known to be NP-complete. Several
heuristic algorithms (such as Gibbs–Poole–Stockmeyer) run in polynomial time while still
producing near-optimal orderings in practice, but exact methods (like
Caprara–Salazar-González) are at least exponential in time complexity and thus are only
feasible for relatively small matrices.

On the other hand, the *matrix bandwidth recognition problem* entails determining whether
there exists a permutation matrix ``P`` such that the bandwidth of ``PAPᵀ`` is at most some
fixed non-negative integer ``k ∈ ℕ``—an optimal permutation that fully minimizes the
bandwidth of ``A`` is not required. Unlike the NP-hard minimization problem, this is
decidable in ``O(nᵏ)`` time.

Many algorithms for both problems exist in the literature, but implementations in the
open-source ecosystem are scarce, with those that do exist primarily tackling older, less
efficient algorithms. This not only makes it difficult for theoretical researchers to
benchmark and compare new approaches but also precludes the application of more performant
alternatives in real-life industry settings. This package aims to bridge this gap,
presenting a unified interface for matrix bandwidth reduction algorithms in Julia. In
addition to providing optimized implementations of many existing approaches,
*MatrixBandwidth.jl* also allows for easy extensibility should researchers wish to test new
ideas, filling a crucial niche in the current research landscape.

The following algorithms are currently supported:
- **Minimization**
    - *Exact*
        - Del Corso–Manzini ([`Minimization.DelCorsoManzini`](@ref))
        - Del Corso–Manzini with perimeter search
            ([`Minimization.DelCorsoManziniWithPS`](@ref))
        - Caprara–Salazar-González ([`Minimization.CapraraSalazarGonzalez`](@ref))
        - Saxe–Gurari–Sudborough ([`Minimization.SaxeGurariSudborough`](@ref))
        - Brute-force search ([`Minimization.BruteForceSearch`](@ref))
    - *Heuristic*
        - Gibbs–Poole–Stockmeyer ([`Minimization.GibbsPooleStockmeyer`](@ref))
        - Cuthill–McKee ([`Minimization.CuthillMcKee`](@ref))
        - Reverse Cuthill–McKee ([`Minimization.ReverseCuthillMcKee`](@ref))
    - *Metaheuristic*
        - Greedy randomized adaptive search procedure (GRASP) ([`Minimization.GRASP`](@ref))
        - Simulated annealing ([`Minimization.SimulatedAnnealing`](@ref))
        - Genetic algorithm ([`Minimization.GeneticAlgorithm`](@ref))
- **Recognition**
    - Del Corso–Manzini ([`Recognition.DelCorsoManzini`](@ref))
    - Del Corso–Manzini with perimeter search ([`Recognition.DelCorsoManziniWithPS`](@ref))
    - Caprara–Salazar-González ([`Recognition.CapraraSalazarGonzalez`](@ref))
    - Saxe–Gurari–Sudborough ([`Recognition.SaxeGurariSudborough`](@ref))
    - Brute-force search ([`Recognition.BruteForceSearch`](@ref))

This package also exports several additional core functions, including (but not limited to)
[`bandwidth`](@ref) and [`profile`](@ref) to compute the original bandwidth and profile of a
matrix prior to any reordering.

The full documentation is available at
[GitHub Pages](https://luis-varona.github.io/MatrixBandwidth.jl/).

# References
- [Maf14](@cite): L. O. Mafteiu-Scai. *The Bandwidths of a Matrix. A Survey of Algorithms*.
    Annals of West University of Timisoara - Mathematics and Computer Science **52**,
    183–223 (2014). https://doi.org/10.2478/awutm-2014-0019.
"""
module MatrixBandwidth

using DataStructures
using Random
using PrecompileTools: @setup_workload, @compile_workload

#= `enqueue!` and `dequeue!` were deprecated in `DataStructures.jl` v0.19 in favor of
`Base.push!` and `Base.popfirst!`, respectively. To maintain backwards compatibility with
older versions of `DataStructures.jl`, we define these methods here if necessary exactly as
they are in the `DataStructures.jl` v0.19 source code. =#
if pkgversion(DataStructures) < v"0.19"
    function Base.push!(q::Queue, x)
        push!(q.store, x)
        return q
    end

    @doc """
        push!(q::Queue, x)

    Inserts the value `x` to the end of the queue `q`.
    """ -> Base.push!

    Base.popfirst!(s::Queue) = popfirst!(s.store)

    @doc """
        popfirst!(q::Queue)

    Removes an element from the front of the queue `q` and returns it.
    """ -> Base.popfirst!
end

include("utils.jl")
include("types.jl")
include("core.jl")

export
    # Submodules
    Minimization,
    Recognition,

    # Types
    AbstractAlgorithm,
    AbstractResult,

    # Core functions
    bandwidth,
    bandwidth_lower_bound,
    profile,

    # Submodule core functions
    minimize_bandwidth,
    has_bandwidth_k_ordering

"""
    const ALGORITHMS :: Dict{Symbol, Union{Dict{Symbol}, Vector}}

A dictionary indexing the data types of all available algorithms by submodule.

For instance, to access all metaheuristic minimization algorithms, use
`MatrixBandwidth.ALGORITHMS[:Minimization][:Metaheuristic]`. Similarly, to access all
recognition algorithms, use `MatrixBandwidth.ALGORITHMS[:Recognition]`.
"""
const ALGORITHMS = Dict{Symbol,Union{Dict{Symbol},Vector}}()

include("Recognition/Recognition.jl")
include("Minimization/Minimization.jl")

using .Minimization, .Recognition

include("startup.jl")

end
