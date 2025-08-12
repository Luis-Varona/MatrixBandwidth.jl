# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestAqua

Static analysis with *Aqua.jl* for the *MatrixBandwidth.jl* package.

*Aqua.jl* offers a general static analyzer, checking for method ambiguities, undefined
`export`s, unbound type parameters, stale dependencies, type piracies, precompilation
issues, and more.
"""
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
