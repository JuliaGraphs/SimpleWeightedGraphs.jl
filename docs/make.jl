using Graphs
using Documenter
using SimpleWeightedGraphs

DocMeta.setdocmeta!(
    SimpleWeightedGraphs, :DocTestSetup, :(using SimpleWeightedGraphs); recursive=true
)

makedocs(;
    modules=[SimpleWeightedGraphs],
    authors="Seth Bromberger and contributors",
    sitename="SimpleWeightedGraphs.jl",
    format=Documenter.HTML(),
    pages=["Home" => "index.md", "Tutorial" => "tutorial.md", "API reference" => "api.md"],
)

deploydocs(; repo="github.com/JuliaGraphs/SimpleWeightedGraphs.jl", devbranch="master")
