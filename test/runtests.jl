# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

using MatrixBandwidth
using Test

const STATIC_ANALYZERS = ["Aqua", "JET"]
const TEST_GROUPS = ["core", "utils"]
const NESTED_TEST_SUITES = ["Minimization/Minimization.jl", "Recognition/Recognition.jl"]

if haskey(ENV, "CI")
    Base.Experimental.suppress_depwarnings(true)
end

# Run static analysis
for analyzer in STATIC_ANALYZERS
    @info "Running static analysis with $analyzer"
    include(joinpath("static_analysis/", "$(lowercase(analyzer)).jl"))
    println()
end

# Run unit tests
for group in TEST_GROUPS
    @info "Testing `$group`"
    include("$group.jl")
    println()
end

# Run more unit tests
for suite in NESTED_TEST_SUITES
    include(suite)
end

# Check that all public names in the package have docstrings
@testset "Docstrings" begin
    if VERSION >= v"1.11" # `Docs.undocumented_names` was only introduced in Julia 1.11
        @info "Checking for undocumented names"
        @test isempty(Docs.undocumented_names(MatrixBandwidth))
    else
        @info "Skipping `Docs.undocumented_names` test: not available on Julia $(VERSION)"
    end
end
