using SearchingTales
using Documenter

DocMeta.setdocmeta!(SearchingTales, :DocTestSetup, :(using SearchingTales); recursive=true)

makedocs(;
    modules=[SearchingTales],
    authors="Laura Brustenga i Moncusí <brust@math.ku.dk> and contributors",
    sitename="SearchingTales.jl",
    format=Documenter.HTML(;
        canonical="https://LauraBMo.github.io/SearchingTales.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LauraBMo/SearchingTales.jl",
    devbranch="main",
)
