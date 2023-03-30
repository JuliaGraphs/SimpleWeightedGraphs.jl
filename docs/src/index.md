```@meta
CurrentModule = SimpleWeightedGraphs
```

# SimpleWeightedGraphs

Documentation for [SimpleWeightedGraphs](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl).

## Quick start

This package defines two new graph types: [`SimpleWeightedGraph`](@ref) and [`SimpleWeightedDiGraph`](@ref).
See the tutorial to discover what you can do with them.
Also refer to the [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) package for more complex algorithms.

## Caveats

Because `SimpleWeighted(Di)Graph`s are stored in sparse matrices, they have two major flaws:

- Iteratively adding or removing vertices or edges is not very efficient. Building the graph in one go from a list of edge sources, destinations and weights is much faster.

- Zero-weight edges are discarded by `add_edge!`. A possible workaround is to [set a very small weight instead](https://stackoverflow.com/questions/48977068/how-to-add-free-edge-to-graph-in-lightgraphs-julia/48994712#48994712).

## Alternatives

If your graphs have more than just edge weights to store, take a look at [MetaGraphsNext.jl](https://github.com/JuliaGraphs/MetaGraphsNext.jl) or [MetaGraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl) for more complex formats.
