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
polynomial-time heuristic approaches (e.g., reverse Cuthill–McKee). Such methods, therefore,
are not feasible for large matrices, but they remain useful when precise solutions are
required for small-to-medium-sized inputs (say, up to ``100×100``).

The following exact algorithms are currently supported:
- Caprara–Salazar-González algorithm ([`CapraraSalazarGonzalez`](@ref))
- Del Corso–Manzini algorithm ([`DelCorsoManzini`](@ref))
- Del Corso–Manzini algorithm with perimeter search ([`DelCorsoManziniWithPS`](@ref))
- Saxe–Gurari–Sudborough algorithm ([`SaxeGurariSudborough`](@ref))

This submodule is part of the `MatrixBandwidth.Minimization` submodule of the
[MatrixBandwidth.jl](https://github.com/Luis-Varona/MatrixBandwidth.jl) package.
"""
module Exact

#! format: off
import ..AbstractSolver, ..NotImplementedError
import .._approach, .._assert_matrix_is_square, .._bool_minimal_band_ordering, .._symmetrize
#! format: on

export CapraraSalazarGonzalez, DelCorsoManzini, DelCorsoManziniWithPS, SaxeGurariSudborough

include("types.jl")

include("solvers/caprara_salazar_gonzalez.jl")
include("solvers/del_corso_manzini.jl")
include("solvers/del_corso_manzini_with_ps.jl")
include("solvers/saxe_gurari_sudborough.jl")

end
