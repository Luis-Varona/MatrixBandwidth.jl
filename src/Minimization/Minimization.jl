# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Minimization

Exact, heuristic, and metaheuristic algorithms for matrix bandwidth minimization in Julia.

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

The following algorithms are currently supported:
- *Exact*
    - Caprara–Salazar-González algorithm ([`CapraraSalazarGonzalez`](@ref))
    - Del Corso–Manzini algorithm ([`DelCorsoManzini`](@ref))
    - Del Corso–Manzini algorithm with perimeter search ([`DelCorsoManziniWithPS`](@ref))
    - Saxe–Gurari–Sudborough algorithm ([`SaxeGurariSudborough`](@ref))
- *Heuristic*
    - Gibbs–Poole–Stockmeyer algorithm ([`GibbsPooleStockmeyer`](@ref))
    - Cuthill–McKee algorithm ([`CuthillMcKee`](@ref))
    - Reverse Cuthill–McKee algorithm ([`ReverseCuthillMcKee`](@ref))
- *Metaheuristic*
    - Greedy randomized adaptive search procedure (GRASP) ([`GRASP`](@ref))
    - Simulated annealing ([`SimulatedAnnealing`](@ref))
    - Genetic algorithm ([`GeneticAlgorithm`](@ref))

This submodule is part of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Minimization

#! format: off
import ..NotImplementedError
import ..bandwidth
import .._assert_matrix_is_square, .._symmetrize
#! format: on

include("types.jl")
include("core.jl")

include("Exact/Exact.jl")
include("Heuristic/Heuristic.jl")
include("Metaheuristic/Metaheuristic.jl")

using .Exact, .Heuristic, .Metaheuristic

# The output struct and core minimization function
export BandMinResult, minimize_bandwidth
export CapraraSalazarGonzalez, # Exact solvers
    DelCorsoManzini,
    DelCorsoManziniWithPS,
    SaxeGurariSudborough
export CuthillMcKee, ReverseCuthillMcKee, GibbsPooleStockmeyer # Heuristic solvers
export SimulatedAnnealing, GeneticAlgorithm, GRASP # Metaheuristic solvers

const DEFAULT_SOLVER = ReverseCuthillMcKee()

end
