# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Recognition

Algorithms for matrix bandwidth recognition in Julia.

The *bandwidth* of an ``n×n`` matrix ``A`` is the minimum non-negative integer
``k ∈ \\{0, 1, …, n - 1\\}`` such that ``A[i, j] = 0`` whenever ``|i - j| > k``. The *matrix
bandwidth recognition problem* entails determining whether there exists a permutation matrix
``P`` such that the bandwidth of ``PAPᵀ`` is at most some fixed non-negative integer (an
optimal permutation that fully minimizes the bandwidth of ``A`` is not required).

The following matrix bandwidth recognition algorithms are currently available:
- Caprara–Salazar-González ([`CapraraSalazarGonzalez`](@ref))
- Del Corso–Manzini ([`DelCorsoManzini`](@ref))
- Del Corso–Manzini with perimeter search ([`DelCorsoManziniWithPS`](@ref))
- Saxe–Gurari–Sudborough ([`SaxeGurariSudborough`](@ref))
- Brute-force search ([`BruteForceSearch`](@ref))

This submodule is part of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Recognition

using MatrixBandwidth
using MatrixBandwidth: NotImplementedError, RectangularMatrixError, StructuralAsymmetryError
using MatrixBandwidth:
    connected_components,
    floyd_warshall_shortest_paths,
    is_structurally_symmetric,
    offdiag_nz_support
using MatrixBandwidth: _requires_structural_symmetry, _problem

using Combinatorics: combinations, permutations
using DataStructures: Queue

export
    # Types
    AbstractDecider,
    RecognitionResult,

    # Core functions
    has_bandwidth_k_ordering,

    # Deciders
    CapraraSalazarGonzalez,
    DelCorsoManzini,
    DelCorsoManziniWithPS,
    SaxeGurariSudborough,
    BruteForceSearch

MatrixBandwidth.ALGORITHMS[:Recognition] = []

include("types.jl")
include("core.jl")

include("deciders/caprara_salazar_gonzalez.jl")
# Defines both `DelCorsoManzini` and `DelCorsoManziniWithPS`
include("deciders/del_corso_manzini.jl")
include("deciders/saxe_gurari_sudborough.jl")
include("deciders/brute_force_search.jl")

const DEFAULT_DECIDER = CapraraSalazarGonzalez()

end
