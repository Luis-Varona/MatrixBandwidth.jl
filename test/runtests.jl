# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

using MatrixBandwidth
using Test

# Run static analysis
for analyzer in readlines(joinpath(@__DIR__, "staticanalyzers"))
    if !isempty(analyzer)
        @info "Running static analysis with $analyzer.jl"
        include("$analyzer.jl")
        println()
    end
end

# Run unit tests
for group in readlines(joinpath(@__DIR__, "testgroups"))
    if !isempty(group)
        @info "Testing `$group`"
        include("$group.jl")
        println()
    end
end

# Check that all public names in the package are documented
@testset "Docstrings" begin
    if VERSION >= v"1.11" # `Docs.undocumented_names` was introduced in Julia 1.11
        @test isempty(Docs.undocumented_names(MatrixBandwidth))
    else
        @info "Skipping `Docs.undocumented_names` test: not available on Julia $(VERSION)"
    end
end
