using CherenkovMediumBase
using Documenter

DocMeta.setdocmeta!(CherenkovMediumBase, :DocTestSetup, :(using CherenkovMediumBase); recursive=true)

makedocs(;
    modules=[CherenkovMediumBase],
    authors="Christian Haack <chr.hck@gmail.com>",
    sitename="CherenkovMediumBase.jl",
    format=Documenter.HTML(;
        canonical="https://juliahep.github.io/CherenkovMediumBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaHEP/CherenkovMediumBase.jl",
    devbranch="main",
)
