using Aqua
using Documenter
using Graphs
using JuliaFormatter
using LinearAlgebra
using SimpleWeightedGraphs
using Test

DocMeta.setdocmeta!(
    SimpleWeightedGraphs, :DocTestSetup, :(using SimpleWeightedGraphs); recursive=true
)

testdir = dirname(@__FILE__)

function testgraphs(g)
    return [g, SimpleWeightedGraph{UInt8,Float64}(g), SimpleWeightedGraph{Int16,Float32}(g)]
end
function testdigraphs(g)
    return [
        g, SimpleWeightedDiGraph{UInt8,Float64}(g), SimpleWeightedDiGraph{Int16,Float32}(g)
    ]
end

testsimplegraphs(g) = [g, Graphs.SimpleGraph{UInt8}(g), Graphs.SimpleGraph{Int16}(g)]
testsimpledigraphs(g) = [g, Graphs.SimpleDiGraph{UInt8}(g), Graphs.SimpleDiGraph{Int16}(g)]

tests = [
    "simpleweightededge",
    "simpleweightedgraph",
    "overrides",
    "persistence",
    "connectivity",
    "a_star",
]

@testset verbose = true "SimpleWeightedGraphs" begin
    @testset verbose = true "Code quality (Aqua.jl)" begin
        Aqua.test_all(SimpleWeightedGraphs; ambiguities=false)
    end
    @testset verbose = false "Code formatting (JuliaFormatter.jl)" begin
        @test format(SimpleWeightedGraphs; verbose=false, overwrite=false)
    end
    @testset verbose = false "Doctests (Documenter.jl)" begin
        doctest(SimpleWeightedGraphs)
    end

    for t in tests
        tp = joinpath(testdir, "$(t).jl")
        include(tp)
    end
end
