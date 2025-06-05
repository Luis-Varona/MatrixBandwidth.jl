# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

module TestAqua

using MatrixBandwidth
using Test
using Aqua

@testset "Static analysis with Aqua" begin
    @test Test.detect_ambiguities(MatrixBandwidth) == Tuple{Method,Method}[]
    Aqua.test_all(MatrixBandwidth)
    @test Aqua.Piracy.hunt(MatrixBandwidth) == Method[]
end

end
