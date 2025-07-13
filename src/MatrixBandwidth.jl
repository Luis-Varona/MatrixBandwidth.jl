# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth

Fast algorithms for matrix bandwidth minimization and matrix bandwidth recognition in Julia.

The *bandwidth* of an ``n×n`` matrix ``A`` is the minimum non-negative integer
``k ∈ [0, n - 1]`` such that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A``
has bandwidth *at most* ``k`` if all entries above the ``k``ᵗʰ superdiagonal and below the
``k``ᵗʰ subdiagonal are zero, and ``A`` has bandwidth *at least* ``k`` if there exists any
nonzero entry in the ``k``ᵗʰ superdiagonal or subdiagonal.

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

The following algorithms are currently supported:
- **Minimization**
    - *Exact*
        - Caprara–Salazar-González algorithm ([`Minimization.CapraraSalazarGonzalez`](@ref))
        - Del Corso–Manzini algorithm ([`Minimization.DelCorsoManzini`](@ref))
        - Del Corso–Manzini algorithm with perimeter search
            ([`Minimization.DelCorsoManziniWithPS`](@ref))
        - Saxe–Gurari–Sudborough algorithm ([`Minimization.SaxeGurariSudborough`](@ref))
        - Brute-force search ([`Minimization.BruteForce`](@ref))
    - *Heuristic*
        - Gibbs–Poole–Stockmeyer algorithm ([`Minimization.GibbsPooleStockmeyer`](@ref))
        - Cuthill–McKee algorithm ([`Minimization.CuthillMcKee`](@ref))
        - Reverse Cuthill–McKee algorithm ([`Minimization.ReverseCuthillMcKee`](@ref))
    - *Metaheuristic*
        - Greedy randomized adaptive search procedure (GRASP) ([`Minimization.GRASP`](@ref))
        - Simulated annealing ([`Minimization.SimulatedAnnealing`](@ref))
        - Genetic algorithm ([`Minimization.GeneticAlgorithm`](@ref))
- **Recognition**
    - Caprara–Salazar-González algorithm ([`Recognition.CapraraSalazarGonzalez`](@ref))
    - Del Corso–Manzini algorithm ([`Recognition.DelCorsoManzini`](@ref))
    - Del Corso–Manzini algorithm with perimeter search
        ([`Recognition.DelCorsoManziniWithPS`](@ref))
    - Saxe–Gurari–Sudborough algorithm ([`Recognition.SaxeGurariSudborough`](@ref))
    - Brute-force search ([`Recognition.BruteForce`](@ref))

[Full documentation](https://Luis-Varona.github.io/MatrixBandwidth.jl/dev/) is available for
the latest development version of this package.
"""
module MatrixBandwidth

using Random

include("utils.jl")
include("types.jl")
include("core.jl")

include("Recognition/Recognition.jl")
include("Minimization/Minimization.jl")

using .Minimization, .Recognition

#= Module exports: allows users to call solvers like `Minimization.GibbsPooleStockmeyer` and
deciders `like Recognition.CapraraSalazarGonzalez`. Solvers/deciders are not exported at the
top level due to name conflicts between `Minimization` and `Recognition`. =#
export Minimization, Recognition

#= Core exports: the original bandwidth (before any reordering) and an `O(n³)` lower bound
from Caprara and Salazar-González (2005). (This bound is tight in many non-trivial cases
but not universally so.) =#
export bandwidth, bandwidth_lower_bound

#= `Minimization` and `Recognition` exports: the core bandwidth minimization and recognition
functions. =#
export minimize_bandwidth, has_bandwidth_k_ordering

# Utility exports: just a random banded matrix generator for now. Useful for test data.
export random_banded_matrix

end
