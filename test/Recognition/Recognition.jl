# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestRecognition

Test suite for the `Recognition` submodule of the `MatrixBandwidth.jl` package.
"""
module TestRecognition

include("core.jl")

include("deciders/caprara_salazar_gonzalez.jl")
include("deciders/del_corso_manzini.jl")
include("deciders/saxe_gurari_sudborough.jl")

end
