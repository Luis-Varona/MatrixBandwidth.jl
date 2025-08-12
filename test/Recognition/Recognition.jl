# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestRecognition

Test suite for the `Recognition` submodule of the *MatrixBandwidth.jl* package.
"""
module TestRecognition

const TEST_GROUPS = [
    "core",
    "deciders/caprara_salazar_gonzalez",
    "deciders/del_corso_manzini",
    "deciders/saxe_gurari_sudborough",
]

for group in TEST_GROUPS
    @info "Testing `Recognition/$group`"
    include("$group.jl")
    println()
end

end
