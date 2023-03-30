```@meta
CurrentModule = SimpleWeightedGraphs
```

# Tutorial

```jldoctest tuto
julia> using Graphs, SimpleWeightedGraphs
```

Here's how to construct an undirected graph (use `SimpleWeightedDiGraph` for directed graphs):

```jldoctest tuto
julia> g = SimpleWeightedGraph(3)
{3, 0} undirected simple Int64 graph with Float64 weights

julia> add_edge!(g, 1, 2, 0.5);

julia> add_edge!(g, 2, 3, 0.8);

julia> add_edge!(g, 1, 3, 2.0);
```

Get the weight of edge from vertex 1 to vertex 2:

```jldoctest tuto
julia> get_weight(g, 1, 2)
0.5
```

Find the shortest path from vertex 1 to vertex 3, taking weights into account:

```jldoctest tuto
julia> enumerate_paths(dijkstra_shortest_paths(g, 1), 3)
3-element Vector{Int64}:
 1
 2
 3
```

Reweight the edge from 1 to 2:

```jldoctest tuto
julia> add_edge!(g, 1, 2, 1.6);
```

Rerun the shortest path calculation from 1 to 3:

```jldoctest tuto
julia> enumerate_paths(dijkstra_shortest_paths(g, 1), 3)
2-element Vector{Int64}:
 1
 3
```

It's possible (and faster) to build the graph from arrays of sources, destinations and weights:

```jldoctest tuto
julia> sources = [1, 2, 1];

julia> destinations = [2, 3, 3];

julia> weights = [0.5, 0.8, 2.0];

julia> g = SimpleWeightedGraph(sources, destinations, weights)
{3, 3} undirected simple Int64 graph with Float64 weights
```

The combine keyword handles repeated pairs (sum by default).
Unexpected results might occur with non-associative combine functions.

```jldoctest tuto
julia> g = SimpleWeightedGraph([1,2,1], [2,1,2], [1,1,1]; combine = +)
{2, 1} undirected simple Int64 graph with Int64 weights

julia> g.weights[2,1] == g.weights[1,2] == 3
true
```

Notice that weights are indexed by `[destination, source]` internally:

```jldoctest tuto
julia> s = SimpleWeightedDiGraph([1,2,1], [2,1,2], [1,1,1]; combine = +);

julia> s.weights[1,2] == 1
true

julia> s.weights[2,1] == 2
true
```
