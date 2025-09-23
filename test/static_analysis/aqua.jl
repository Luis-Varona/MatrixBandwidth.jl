# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestAqua

Static analysis with Aqua.jl for the MatrixBandwidth.jl package.

Aqua.jl offers a general static analyzer, checking for method ambiguities, undefined
`export`s, unbound type parameters, stale dependencies, type piracies, precompilation
issues, and more.
"""
module TestAqua

using MatrixBandwidth
using Test
using Aqua

@testset "Static analysis with Aqua" begin
    @test Test.detect_ambiguities(MatrixBandwidth) == Tuple{Method,Method}[]
    Aqua.test_all(
        MatrixBandwidth,
        piracies=(; treat_as_own=[Base.push!, Base.popfirst!]),
        # Account for our manual definitions of `Base.push!` and `Base.popfirst!`
        persistent_tasks=false,
    )
    #= We manually define the `Base.push!(::Queue, ::Any)` and `Base.popfirst!(::Queue)`
    methods when using `DataStructures.jl` v0.18, as `enqueue!` and `dequeue!` were only
    deprecated in v0.19+ and we need to maintain compatibility with v0.18 as well. =#
    @test length(Aqua.Piracy.hunt(MatrixBandwidth)) <= 2
end

end
