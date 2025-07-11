# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestMinimization

Test suite for the `Minimization` submodule of the `MatrixBandwidth.jl` package.
"""
module TestMinimization

include("core.jl")

include("Exact/Exact.jl")
include("Heuristic/Heuristic.jl")
include("Metaheuristic/Metaheuristic.jl")

end
