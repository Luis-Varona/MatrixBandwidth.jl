# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Minimization.Exact

Exact solvers for matrix bandwidth minimization.

Exact methods are those which guarantee an optimal ordering producing the true minimum
bandwidth of a matrix. Since bandwidth minimization is an NP-complete problem, existing
exact algorithms are, at best, exponential in time complexity—much worse than many
polynomial-time heuristic approaches (e.g., Gibbs–Poole–Stockmeyer). Such methods,
therefore, are not feasible for large matrices, but they remain useful when precise
solutions are required for small-to-medium-sized inputs (say, up to ``100×100``).

The following exact matrix bandwidth minimization algorithms are currently available:
- Del Corso–Manzini ([`DelCorsoManzini`](@ref))
- Del Corso–Manzini with perimeter search ([`DelCorsoManziniWithPS`](@ref))
- Caprara–Salazar-González ([`CapraraSalazarGonzalez`](@ref))
- Saxe–Gurari–Sudborough ([`SaxeGurariSudborough`](@ref))
- Brute-force search ([`BruteForceSearch`](@ref))

This submodule is part of the `MatrixBandwidth.Minimization` submodule of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Exact

using MatrixBandwidth
using MatrixBandwidth: NotImplementedError, StructuralAsymmetryError
using MatrixBandwidth: connected_components, floyd_warshall_shortest_paths
using MatrixBandwidth: _requires_structural_symmetry

using MatrixBandwidth.Recognition

using MatrixBandwidth.Minimization
using MatrixBandwidth.Minimization: _approach, _minimize_bandwidth_impl

using Combinatorics

export
    # Types
    ExactSolver,

    # Exact solvers
    DelCorsoManzini,
    DelCorsoManziniWithPS,
    CapraraSalazarGonzalez,
    SaxeGurariSudborough,
    BruteForceSearch

MatrixBandwidth.ALGORITHMS[:Minimization][:Exact] = []

include("types.jl")

# Defines both `DelCorsoManzini` and `DelCorsoManziniWithPS`
include("solvers/del_corso_manzini.jl")
include("solvers/caprara_salazar_gonzalez.jl")
include("solvers/saxe_gurari_sudborough.jl")
include("solvers/brute_force_search.jl")

end
