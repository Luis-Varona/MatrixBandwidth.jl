# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

module TestJET

using MatrixBandwidth
using Test
using JET

@testset "Static analysis with JET" begin
    rep = report_package("MatrixBandwidth")
    jet_reports = JET.get_reports(rep)

    @show length(jet_reports)
    @show rep

    @test length(jet_reports) < 30
    @test_broken length(jet_reports) == 0
end

end
