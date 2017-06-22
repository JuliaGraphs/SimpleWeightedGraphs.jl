# SimpleWeightedGraphs

[![Build Status](https://travis-ci.org/sbromberger/SimpleWeightedGraphs.jl.svg?branch=master)](https://travis-ci.org/sbromberger/SimpleWeightedGraphs.jl)
[![codecov.io](http://codecov.io/github/sbromberger/SimpleWeightedGraphs.jl/coverage.svg?branch=master)](http://codecov.io/github/sbromberger/SimpleWeightedGraphs.jl?branch=master)

Weighted Graphs for [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl).

Usage:
```
using LightGraphs, SimpleWeightedGraphs

g = SimpleWeightedGraph(3)
add_edge!(g, 1, 2, 0.5)
add_edge!(g, 2, 3, 0.8)
add_edge!(g, 1, 3, 2.0)

# find the shortest path from vertex 1 to vertex 3 taking weights into account.
enumerate_paths(dijkstra_shortest_paths(g, 1), 3)
3-element Array{Int64,1}:
 1
 2
 3
```
