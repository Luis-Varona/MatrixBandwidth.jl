using MatrixBandwidth
using Documenter

DocMeta.setdocmeta!(
    MatrixBandwidth, :DocTestSetup, :(using MatrixBandwidth); recursive=true
)

makedocs(;
    modules=[MatrixBandwidth],
    authors="Luis M. B. Varona <lbvarona@mta.ca>",
    sitename="MatrixBandwidth.jl",
    format=Documenter.HTML(;
        canonical="https://Luis-Varona.github.io/MatrixBandwidth.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/Luis-Varona/MatrixBandwidth.jl", devbranch="main")
