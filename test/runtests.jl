# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

using MatrixBandwidth
using Test

for analyzer in readlines(joinpath(@__DIR__, "staticanalyzers"))
    @info "Running static analysis with $analyzer"
    include("static_analysis/$(lowercase(analyzer)).jl")
    println()
end

for file in readlines(joinpath(@__DIR__, "testgroups"))
    @info "Testing $file"
    include("$file.jl")
    println()
end

@testset "Docstrings" begin
    if VERSION >= v"1.11" # `Docs.undocumented_names` was introduced in Julia 1.11
        @test isempty(Docs.undocumented_names(MatrixBandwidth))
    else
        @info "Skipping `Docs.undocumented_names` test: not available on Julia $(VERSION)"
    end
end
