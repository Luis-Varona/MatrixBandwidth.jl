# Copyright 2025 Luis M. B. Varona
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

using MatrixBandwidth
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(
    MatrixBandwidth, :DocTestSetup, :(using MatrixBandwidth); recursive=true
)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:alpha)

makedocs(;
    modules=[MatrixBandwidth, MatrixBandwidth.Minimization, MatrixBandwidth.Recognition],
    authors="Luis M. B. Varona <lbvarona@mta.ca>",
    sitename="MatrixBandwidth.jl",
    format=Documenter.HTML(;
        canonical="https://Luis-Varona.github.io/MatrixBandwidth.jl",
        edit_link="main",
        assets=["assets/styles.css"],
        size_threshold=1_000_000,
        size_threshold_warn=200_000,
    ),
    plugins=[bib],
    pages=[
        "Home" => "index.md",
        "Public API" => "public_api.md",
        "Private API" => "private_api.md",
    ],
)

deploydocs(; repo="github.com/Luis-Varona/MatrixBandwidth.jl", devbranch="main")
