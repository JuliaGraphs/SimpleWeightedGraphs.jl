using Graphs
using SimpleWeightedGraphs
using LinearAlgebra
using Test

testdir = dirname(@__FILE__)

testgraphs(g) = [g, SimpleWeightedGraph{UInt8,Float64}(g), SimpleWeightedGraph{Int16,Float32}(g)]
testdigraphs(g) = [g, SimpleWeightedDiGraph{UInt8,Float64}(g), SimpleWeightedDiGraph{Int16,Float32}(g)]

testsimplegraphs(g) = [g, Graphs.SimpleGraph{UInt8}(g), Graphs.SimpleGraph{Int16}(g)]
testsimpledigraphs(g) = [g, Graphs.SimpleDiGraph{UInt8}(g), Graphs.SimpleDiGraph{Int16}(g)]

tests = [
    "simpleweightededge",
    "simpleweightedgraph",
    "overrides",
    "persistence",
    "connectivity",
    "a_star"
]

@testset "SimpleWeightedGraphs" begin
    for t in tests
        tp = joinpath(testdir, "$(t).jl")
        include(tp)
    end
end
