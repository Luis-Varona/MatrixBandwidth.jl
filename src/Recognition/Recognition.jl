# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Recognition

Algorithms for matrix bandwidth recognition in Julia.

The *bandwidth* of an ``n×n`` matrix ``A`` is the minimum non-negative integer
``k ∈ \\{0, 1, …, n - 1\\}`` such that ``A[i, j] = 0`` whenever ``|i - j| > k``.
Equivalently, ``A`` has bandwidth *at most* ``k`` if all entries above the ``k``ᵗʰ
superdiagonal and below the ``k``ᵗʰ subdiagonal are zero, and ``A`` has bandwidth *at least*
``k`` if there exists any nonzero entry in the ``k``ᵗʰ superdiagonal or subdiagonal.

The *matrix bandwidth recognition problem* entails determining whether there exists a
permutation matrix ``P`` such that the bandwidth of ``PAPᵀ`` is at most some fixed
non-negative integer ``k ∈ ℕ``—an optimal permutation that fully minimizes the bandwidth of
``A`` is not required. Unlike the NP-hard minimization problem, this is decidable in
``O(nᵏ)`` time.

The following algorithms are currently supported:
- Caprara–Salazar-González algorithm ([`CapraraSalazarGonzalez`](@ref))
- Del Corso–Manzini algorithm ([`DelCorsoManzini`](@ref))
- Del Corso–Manzini algorithm with perimeter search ([`DelCorsoManziniWithPS`](@ref))
- Saxe–Gurari–Sudborough algorithm ([`SaxeGurariSudborough`](@ref))
- Brute-force search ([`BruteForceSearch`](@ref))

This submodule is part of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Recognition

#! format: off
import ..ALGORITHMS
import ..AbstractAlgorithm, ..AbstractResult
import ..NotImplementedError, ..RectangularMatrixError, ..StructuralAsymmetryError
import ..bandwidth, ..bandwidth_lower_bound
import .._problem
import .._connected_components, .._is_structurally_symmetric, .._offdiag_nonzero_support
#! format: on

using Combinatorics: combinations, permutations
using DataStructures: Queue

# THe output struct and core recognition function
export RecognitionResult, has_bandwidth_k_ordering
export CapraraSalazarGonzalez, # Recognition algorithms
    DelCorsoManzini,
    DelCorsoManziniWithPS,
    SaxeGurariSudborough,
    BruteForceSearch

ALGORITHMS[:Recognition] = []

include("types.jl")
include("core.jl")

include("deciders/caprara_salazar_gonzalez.jl")
# Defines both `DelCorsoManzini` and `DelCorsoManziniWithPS`
include("deciders/del_corso_manzini.jl")
include("deciders/saxe_gurari_sudborough.jl")
include("deciders/brute_force_search.jl")

const DEFAULT_DECIDER = CapraraSalazarGonzalez()

end
