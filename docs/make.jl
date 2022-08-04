using GeneralizedMultiparticleMie
using Documenter

DocMeta.setdocmeta!(GeneralizedMultiparticleMie, :DocTestSetup, :(using GeneralizedMultiparticleMie); recursive=true)

makedocs(;
    modules=[GeneralizedMultiparticleMie],
    authors="Gabriel Wu <wuzihua@pku.edu.cn> and contributors",
    repo="https://github.com/lucifer1004/GeneralizedMultiparticleMie.jl/blob/{commit}{path}#{line}",
    sitename="GeneralizedMultiparticleMie.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lucifer1004.github.io/GeneralizedMultiparticleMie.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lucifer1004/GeneralizedMultiparticleMie.jl",
    devbranch="main",
)
