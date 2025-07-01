# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    MatrixBandwidth.Recognition

Algorithms for matrix bandwidth recognition in Julia.

The *bandwidth* of a square matrix ``A`` is the minimum non-negative integer ``k ∈ ℕ`` such
that ``A[i, j] = 0`` whenever ``|i - j| > k``. Equivalently, ``A`` has bandwidth *at most*
``k`` if all entries above the ``k``-th superdiagonal and below the ``k``-th subdiagonal are
zero, and ``A`` has bandwidth *at least* ``k`` if there exists any nonzero entry in the
``k``-th superdiagonal or subdiagonal.

The *matrix bandwidth recognition problem* [TODO: Write here]

The following algorithms are currently supported [TODO: Add refs later]:
- Caprara–Salazar-González algorithm
- Del Corso–Manzini algorithm
- Del Corso–Manzini algorithm with perimeter search
- Saxe–Gurari–Sudborough algorithm

This submodule is part of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Recognition
# TODO: Add comments

#! format: off
import ..NotImplementedError
import ..bandwidth
import .._symmetrize
#! format: on

# include("types.jl")
# include("core.jl")

# include("deciders/caprara_salazar_gonzalez.jl")
# include("deciders/del_corso_manzini.jl")
# include("deciders/del_corso_manzini_with_ps.jl")
# include("deciders/saxe_gurari_sudborough.jl")

# export BandRecogResult, has_bandwidth_k_ordering
# export CapraraSalazarGonzalez, # Recognition algorithms
#     DelCorsoManzini,
#     DelCorsoManziniWithPS,
#     SaxeGurariSudborough

# const DEFAULT_DECIDER = CapraraSalazarGonzalez()

end
