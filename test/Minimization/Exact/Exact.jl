# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestExact

Test suite for the `Minimization.Exact` submodule of the `MatrixBandwidth.jl` package.
"""
module TestExact

include("solvers/caprara_salazar_gonzalez.jl")
include("solvers/del_corso_manzini.jl")
include("solvers/saxe_gurari_sudborough.jl")

end
