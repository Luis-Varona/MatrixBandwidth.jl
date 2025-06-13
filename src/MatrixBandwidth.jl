# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth

Exact, heuristic, and metaheuristic algorithms for matrix bandwidth minimization in Julia.

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

*Matrix bandwidth minimization* is the problem of finding a permutation matrix ``P`` so that
the bandwidth of ``PAPᵀ`` is minimized; this is known to be NP-complete. Several heuristic
algorithms (such as reverse Cuthill–McKee) run in polynomial time while still producing
near-optimal orderings in practice, but exact methods (like MB-PS) are exponential in time
complexity and thus only feasible for relatively small matrices.

The following algorithms are currently supported:
- **Exact**
    - [`MBID`](@ref): Minimum bandwidth by iterative deepening (MB-ID)
    - [`MBPS`](@ref): Minimum bandwidth by perimeter search (MB-PS)
- **Heuristic**
    - [`CuthillMcKee`](@ref): Cuthill–McKee algorithm
    - [`ReverseCuthillMcKee`](@ref): Reverse Cuthill–McKee algorithm
- **Metaheuristic**
    - [`SimulatedAnnealing`](@ref): Simulated annealing
    - [`GeneticAlgorithm`](@ref): Genetic algorithm
    - [`GRASP`](@ref): Greedy randomized adaptive search procedure (GRASP)

[Full documentation](https://Luis-Varona.github.io/MatrixBandwidth.jl/dev/) is available for
the latest development version of this package.
"""
module MatrixBandwidth

using Random

include("utils.jl")
include("types.jl")
include("core.jl")

include("Exact/Exact.jl")
include("Heuristic/Heuristic.jl")
include("Metaheuristic/Metaheuristic.jl")

using .Exact, .Heuristic, .Metaheuristic

# The output struct, main minimization function, and raw matrix bandwidth function
export BandwidthResult, minimize_bandwidth, bandwidth
export MBID, MBPS # Exact solvers
export CuthillMcKee, ReverseCuthillMcKee # Heuristic solvers
export SimulatedAnnealing, GeneticAlgorithm, GRASP # Metaheuristic solvers
export random_banded_matrix # Random banded matrices for test data

const DEFAULT_SOLVER = ReverseCuthillMcKee()

end
