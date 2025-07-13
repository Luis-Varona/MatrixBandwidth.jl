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

const TEST_GROUPS = [
    "solvers/caprara_salazar_gonzalez",
    "solvers/del_corso_manzini",
    "solvers/saxe_gurari_sudborough",
]

for group in TEST_GROUPS
    @info "Testing `Minimization/Exact/$group`"
    include("$group.jl")
    println()
end

end
