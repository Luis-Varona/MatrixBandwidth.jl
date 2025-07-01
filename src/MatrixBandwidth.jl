# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth

[TODO: Update to reflect recognition algorithms too]

Fast algorithms for matrix bandwidth minimization and matrix bandwidth recognition in Julia.

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

The *matrix bandwidth minimization problem* entails finding a permutation matrix ``P`` so
that the bandwidth of ``PAPᵀ`` is minimized; this is known to be NP-complete. Several
heuristic algorithms (such as reverse Cuthill–McKee) run in polynomial time while still
producing near-optimal orderings in practice, but exact methods (like
Caprara–Salazar-González) are exponential in time complexity and thus are only feasible for
relatively small matrices.

On the other hand, the *matrix bandwidth recognition problem* [TODO: Write here]

The following algorithms are currently supported [TODO: Add refs for `Recognition` later]:
- **Minimization**
    - *Exact*
        - Caprara–Salazar-González algorithm ([`CapraraSalazarGonzalez`](@ref))
        - Del Corso–Manzini algorithm ([`DelCorsoManzini`](@ref))
        - Del Corso–Manzini algorithm with perimeter search
          ([`DelCorsoManziniWithPS`](@ref))
        - Saxe–Gurari–Sudborough algorithm ([`SaxeGurariSudborough`](@ref))
    - *Heuristic*
        - Gibbs–Poole–Stockmeyer algorithm ([`GibbsPooleStockmeyer`](@ref))
        - Cuthill–McKee algorithm ([`CuthillMcKee`](@ref))
        - Reverse Cuthill–McKee algorithm ([`ReverseCuthillMcKee`](@ref))
    - *Metaheuristic*
        - Greedy randomized adaptive search procedure (GRASP) ([`GRASP`](@ref))
        - Simulated annealing ([`SimulatedAnnealing`](@ref))
        - Genetic algorithm ([`GeneticAlgorithm`](@ref))
- **Recognition**
    - Caprara–Salazar-González algorithm
    - Del Corso–Manzini algorithm
    - Del Corso–Manzini algorithm with perimeter search
    - Saxe–Gurari–Sudborough algorithm

[Full documentation](https://Luis-Varona.github.io/MatrixBandwidth.jl/dev/) is available for
the latest development version of this package.
"""
module MatrixBandwidth

using Random

include("utils.jl")
include("core.jl")

include("Minimization/Minimization.jl")
include("Recognition/Recognition.jl")

using .Minimization, .Recognition

export Minimization, Recognition # TODO: Comment here
export bandwidth # Raw matrix bandwidth computation (with the current permutation)
export BandMinResult, minimize_bandwidth # Matrix bandwidth minimization
# export BandRecogResult, has_bandwidth_k_ordering # Matrix bandwidth recognition
export random_banded_matrix # Random banded matrix generation for test data

end
