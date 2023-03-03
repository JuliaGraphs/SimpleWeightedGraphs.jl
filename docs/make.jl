using SimpleWeightedGraphs
using Documenter

DocMeta.setdocmeta!(SimpleWeightedGraphs, :DocTestSetup, :(using SimpleWeightedGraphs); recursive=true)

makedocs(;
    modules=[SimpleWeightedGraphs],
    authors="Seth Bromberger and contributors",
    repo="https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl/blob/{commit}{path}#{line}",
    sitename="SimpleWeightedGraphs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGraphs.github.io/SimpleWeightedGraphs.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API reference" => "api.md",
    ],
    linkcheck=true,
    strict=true,
)

deploydocs(;
    repo="github.com/JuliaGraphs/SimpleWeightedGraphs.jl",
    devbranch="master",
)
