var documenterSearchIndex = {"docs":
[{"location":"api/#API-reference","page":"API reference","title":"API reference","text":"","category":"section"},{"location":"api/#Index","page":"API reference","title":"Index","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"","category":"page"},{"location":"api/#Docstrings","page":"API reference","title":"Docstrings","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [SimpleWeightedGraphs]","category":"page"},{"location":"api/#SimpleWeightedGraphs.SimpleWeightedGraphs","page":"API reference","title":"SimpleWeightedGraphs.SimpleWeightedGraphs","text":"SimpleWeightedGraphs\n\nA package for graphs with edge weights and no self-loops, stored as sparse adjacency matrices.\n\n\n\n\n\n","category":"module"},{"location":"api/#SimpleWeightedGraphs.AbstractSimpleWeightedEdge","page":"API reference","title":"SimpleWeightedGraphs.AbstractSimpleWeightedEdge","text":"AbstractSimpleWeightedEdge{T}\n\nAbstract type for weighted edges with endpoints of type T.\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleWeightedGraphs.AbstractSimpleWeightedGraph","page":"API reference","title":"SimpleWeightedGraphs.AbstractSimpleWeightedGraph","text":"AbstractSimpleWeightedGraph\n\nAn abstract type representing a simple graph structure with edge weights.\n\nThe only requirement for concrete subtypes is that they should implement a method Graphs.weights(g) which returns the weighted adjacency matrix.\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleWeightedGraphs.SWGFormat","page":"API reference","title":"SimpleWeightedGraphs.SWGFormat","text":"SWGFormat\n\nThe storage format of SimpleWeightedGraph files. Multiple graphs may be present in one file.\n\nFor each graph, the file contains a one line header:\n\n\"LightGraphs.SimpleWeightedGraph\", <num_vertices>, <num_edges>, {\"d\" | \"u\"}, <name>[, <ver>, <vdatatype>, <wdatatype>, <graphcode>]\n\n\"LightGraphs.SimpleWeightedGraph\" is a fixed string\n<num_vertices> is an integer\n<num_edges> is an integer\n\"d\" for directed graph, \"u\" for undirected (note that this   option does not perform any additional edge construction, it's   merely used to return the correct type of graph)\n<name> is a string\n<ver> is an int\n<vdatatype> is a string (\"UInt8\", etc.)\n<wdatatype> is a string describing the data type of the weights\n<graphcode> is a string.\n\nThe header is followed by a list of (comma-delimited) edges, each on a separate line:\n\n<src>,<dst>,<weight>\n\n<src> is an int giving the source of the edge\n<dst> is an int giving the destination of the edge\n<weight> is a real giving the weight of the edge\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleWeightedGraphs.SimpleWeightedDiGraph","page":"API reference","title":"SimpleWeightedGraphs.SimpleWeightedDiGraph","text":"SimpleWeightedDiGraph{T,U}\n\nA type representing a directed weighted graph with vertices of type T and edge weights of type U.\n\nFields\n\nweights::SparseMatrixCSC{U,T}: weighted adjacency matrix, indexed by (dst, src)\n\ntip: Performance\nIteratively adding/removing vertices or edges is not very efficient for this type of graph: better construct the graph in one shot if possible.\n\nBasic constructors\n\nSimpleWeightedDiGraph()  # empty\nSimpleWeightedDiGraph(n)  # n vertices, no edges\nSimpleWeightedDiGraph(graph)  # from graph\nSimpleWeightedDiGraph(adjmx; permute)  # from adjacency matrix, possibly transposed\nSimpleWeightedDiGraph(sources, destinations, weights)  # from list of edges\n\nUse methods(SimpleWeightedDiGraph) for the full list of constructors. When building a new graph from a list of edges, be aware that repeating (src, dst) pairs may lead to undefined behavior (e.g. due to floating point errors during weight addition).\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleWeightedGraphs.SimpleWeightedDiGraphEdge","page":"API reference","title":"SimpleWeightedGraphs.SimpleWeightedDiGraphEdge","text":"SimpleWeightedDiGraphEdge\n\nAlias for SimpleWeightedEdge.\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleWeightedGraphs.SimpleWeightedEdge","page":"API reference","title":"SimpleWeightedGraphs.SimpleWeightedEdge","text":"SimpleWeightedEdge{T,U}\n\nConcrete struct for a weighted edge with endpoints of type T and a weight of type U<:Real.\n\nFields\n\nsrc::T: edge source\ndst::T: edge destination\nweight::U: edge weight\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleWeightedGraphs.SimpleWeightedEdge-Tuple{Any, Any}","page":"API reference","title":"SimpleWeightedGraphs.SimpleWeightedEdge","text":"SimpleWeightedEdge(u, v)\n\nConstruct a SimpleWeightedEdge from u to v with a default weight of 1.0.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.SimpleWeightedEdge-Tuple{Pair}","page":"API reference","title":"SimpleWeightedGraphs.SimpleWeightedEdge","text":"SimpleWeightedEdge(u => v)\n\nConstruct a SimpleWeightedEdge from u to v with a default weight of 1.0.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.SimpleWeightedEdge-Tuple{Tuple{T, T, T} where T}","page":"API reference","title":"SimpleWeightedGraphs.SimpleWeightedEdge","text":"SimpleWeightedEdge((u, v, w))\n\nConstruct a SimpleWeightedEdge from u to v with a weight of w.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.SimpleWeightedEdge-Tuple{Tuple{T, T} where T}","page":"API reference","title":"SimpleWeightedGraphs.SimpleWeightedEdge","text":"SimpleWeightedEdge((u, v))\n\nConstruct a SimpleWeightedEdge from u to v with a default weight of 1.0.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.SimpleWeightedGraph","page":"API reference","title":"SimpleWeightedGraphs.SimpleWeightedGraph","text":"SimpleWeightedGraph{T, U}\n\nA type representing an undirected weighted graph with vertices of type T and edge weights of type U.\n\nFields\n\nweights::SparseMatrixCSC{U,T}: weighted adjacency matrix, indexed by (dst, src)\n\ntip: Performance\nIteratively adding/removing vertices or edges is not very efficient for this type of graph: better construct the graph in one shot if possible.\n\nBasic constructors\n\nSimpleWeightedGraph()  # empty\nSimpleWeightedGraph(n)  # n vertices, no edges\nSimpleWeightedGraph(graph)  # from graph\nSimpleWeightedGraph(adjmx)  # from adjacency matrix\nSimpleWeightedGraph(sources, destinations, weights)  # from list of edges\n\nUse methods(SimpleWeightedGraph) for the full list of constructors. When building a new graph from a list of edges, be aware that repeating (src, dst) pairs may lead to undefined behavior (e.g. due to floating point errors during weight addition).\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleWeightedGraphs.SimpleWeightedGraphEdge","page":"API reference","title":"SimpleWeightedGraphs.SimpleWeightedGraphEdge","text":"SimpleWeightedGraphEdge\n\nAlias for SimpleWeightedEdge.\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleWeightedGraphs.WDiGraph","page":"API reference","title":"SimpleWeightedGraphs.WDiGraph","text":"WDiGraph\n\nAlias for SimpleWeightedDiGraph.\n\n\n\n\n\n","category":"type"},{"location":"api/#SimpleWeightedGraphs.WGraph","page":"API reference","title":"SimpleWeightedGraphs.WGraph","text":"WGraph\n\nAlias for SimpleWeightedGraph.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.getindex-Tuple{SimpleWeightedDiGraph, AbstractEdge, Val{:weight}}","page":"API reference","title":"Base.getindex","text":"g[e, Val(:weight)]\n\nReturn the weight of edge e.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.getindex-Tuple{SimpleWeightedDiGraph, Integer, Integer, Val{:weight}}","page":"API reference","title":"Base.getindex","text":"g[i, j, Val(:weight)]\n\nReturn the weight of edge (i, j).\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.getindex-Tuple{SimpleWeightedGraph, AbstractEdge, Val{:weight}}","page":"API reference","title":"Base.getindex","text":"g[e, :weight]\n\nReturn the weight of edge e.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.getindex-Tuple{SimpleWeightedGraph, Integer, Integer, Val{:weight}}","page":"API reference","title":"Base.getindex","text":"g[i, j, :weight]\n\nReturn the weight of edge (i, j).\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.LinAlg.adjacency_matrix","page":"API reference","title":"Graphs.LinAlg.adjacency_matrix","text":"Graphs.adjacency_matrix(g, T; dir)\n\nConstruct the weighted adjacency matrix, filled with element type T and considering edge direction dir ∈ [:in, :out, :both] (default is :out).\n\n\n\n\n\n","category":"function"},{"location":"api/#Graphs.LinAlg.laplacian_matrix","page":"API reference","title":"Graphs.LinAlg.laplacian_matrix","text":"Graphs.laplacian_matrix(g, T; dir)\n\nSubtract the adjacency matrix to the degree matrix, both filled with element type T and considering edge direction dir ∈ [:in, :out, :both] (unlike in Graphs.jl, default is :out).\n\n\n\n\n\n","category":"function"},{"location":"api/#Graphs.SimpleGraphs.add_edge!-Union{Tuple{U}, Tuple{T}, Tuple{AbstractSimpleWeightedGraph{T, U}, Any, Any, Any}} where {T, U}","page":"API reference","title":"Graphs.SimpleGraphs.add_edge!","text":"Graphs.add_edge!(g, u, v, w)\n\nAdd the edge (u, v) to the graph with a weight of w.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.SimpleGraphs.add_edge!-Union{Tuple{U}, Tuple{T}, Tuple{AbstractSimpleWeightedGraph{T, U}, Any, Any}} where {T, U}","page":"API reference","title":"Graphs.SimpleGraphs.add_edge!","text":"Graphs.add_edge!(g, u, v)\n\nAdd the edge (u, v) to the graph with a default weight of 1.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.SimpleGraphs.add_edge!-Union{Tuple{U}, Tuple{T}, Tuple{AbstractSimpleWeightedGraph{T, U}, Any}} where {T, U}","page":"API reference","title":"Graphs.SimpleGraphs.add_edge!","text":"Graphs.add_edge!(g, e)\n\nAdd the edge e to the graph, where e is any object that can be converted into edgetype(g).\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.SimpleGraphs.rem_edge!-Tuple{SimpleWeightedDiGraph, AbstractEdge}","page":"API reference","title":"Graphs.SimpleGraphs.rem_edge!","text":"Graphs.rem_edge!(g, e)\n\nRemove the edge e from the graph.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.SimpleGraphs.rem_edge!-Tuple{SimpleWeightedGraph, AbstractEdge}","page":"API reference","title":"Graphs.SimpleGraphs.rem_edge!","text":"Graphs.rem_edge!(g, e)\n\nRemove the edge e from the graph.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.SimpleGraphs.rem_edge!-Union{Tuple{U}, Tuple{T}, Tuple{AbstractSimpleWeightedGraph{T, U}, Integer, Integer}} where {T, U}","page":"API reference","title":"Graphs.SimpleGraphs.rem_edge!","text":"Graphs.rem_edge!(g, u, v)\n\nRemove the edge (u, v) from the graph.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.SimpleGraphs.rem_vertex!-Tuple{SimpleWeightedDiGraph, Integer}","page":"API reference","title":"Graphs.SimpleGraphs.rem_vertex!","text":"rem_vertex!(g::SimpleWeightedDiGraph, v)\n\nRemove the vertex v from graph g. Return false if removal fails (e.g., if vertex is not in the graph) and true otherwise.\n\ntip: Correctness\nThis operation has to be performed carefully if one keeps external data structures indexed by edges or vertices in the graph, since internally the removal results in all vertices with indices greater than v being shifted down one.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.SimpleGraphs.rem_vertex!-Tuple{SimpleWeightedGraph, Integer}","page":"API reference","title":"Graphs.SimpleGraphs.rem_vertex!","text":"rem_vertex!(g::SimpleWeightedGraph, v)\n\nRemove the vertex v from graph g. Return false if removal fails (e.g., if vertex is not in the graph) and true otherwise.\n\ntip: Correctness\nThis operation has to be performed carefully if one keeps external data structures indexed by edges or vertices in the graph, since internally the removal results in all vertices with indices greater than v being shifted down one.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.cartesian_product-Union{Tuple{G}, Tuple{G, G}} where G<:AbstractSimpleWeightedGraph","page":"API reference","title":"Graphs.cartesian_product","text":"Graphs.cartesian_product(g, h)\n\nCompute the weighted cartesian product of two weighted graphs.\n\nwarning: Warning\nIt is possible that this is suboptimal, but it is the most trivial extension of the implementation used in Graphs.jl.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.connected_components-Union{Tuple{SimpleWeightedGraph{T, U}}, Tuple{U}, Tuple{T}} where {T, U}","page":"API reference","title":"Graphs.connected_components","text":"Graphs.connected_components(g)\n\nCompute the connected components of a weighted graph. Note that an edge with weight 0 will still be counted as an edge if it exists in the sparse weights matrix.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.has_edge-Union{Tuple{U}, Tuple{T}, Tuple{AbstractSimpleWeightedGraph{T, U}, Any, Any}} where {T, U}","page":"API reference","title":"Graphs.has_edge","text":"Graphs.has_edge(g, u, v)\n\nCheck the existence of the edge (u, v) in the graph.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.has_edge-Union{Tuple{U}, Tuple{T}, Tuple{AbstractSimpleWeightedGraph{T, U}, Any}} where {T, U}","page":"API reference","title":"Graphs.has_edge","text":"Graphs.has_edge(g, e)\n\nCheck the existence of the edge e in the graph, where e is any object that can be converted into edgetype(g).\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.induced_subgraph-Union{Tuple{U}, Tuple{G}, Tuple{G, AbstractVector{U}}} where {G<:AbstractSimpleWeightedGraph, U<:Integer}","page":"API reference","title":"Graphs.induced_subgraph","text":"Graphs.induced_subgraph(g, vlist)\n\nCompute the weighted subgraph induced by a list of vertices.\n\nReturn a tuple containing the new graph and the list of vertices.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.inneighbors-Tuple{SimpleWeightedDiGraph, Integer}","page":"API reference","title":"Graphs.inneighbors","text":"Graphs.inneighbors(g::SimpleWeightedDiGraph, v)\n\nReturn the vector of inneighbors of vertex v.\n\ntip: Performance\nThis function is less efficient than inneighbors for directed weighted graphs (it allocates a new vector).\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.loadgraph-Tuple{IO, String, SWGFormat}","page":"API reference","title":"Graphs.loadgraph","text":"Graphs.loadgraph(io, gname, SWGFormat())\n\nReturn a single graph with name gname loaded from the IO stream io using SWGFormat.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.loadgraphs-Tuple{IO, SWGFormat}","page":"API reference","title":"Graphs.loadgraphs","text":"Graphs.loadgraphs(io, SWGFormat())\n\nReturn a dictionary of name => graph pairs loaded from the IO stream io using SWGFormat.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.outneighbors-Tuple{SimpleWeightedDiGraph, Integer}","page":"API reference","title":"Graphs.outneighbors","text":"Graphs.outneighbors(g::SimpleWeightedDiGraph, v)\n\nReturn the vector of outneighbors of vertex v.\n\ntip: Performance\nThis function is more efficient than inneighbors for directed weighted graphs.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.pagerank","page":"API reference","title":"Graphs.pagerank","text":"Graphs.pagerank(g, α=0.85, n=100, ϵ=1.0e-6)\n\nApply the page rank algorithm on a weighted graph.\n\n\n\n\n\n","category":"function"},{"location":"api/#Graphs.savegraph","page":"API reference","title":"Graphs.savegraph","text":"Graphs.savegraph(fn, g, gname)\n\nwarning: Warning\nThis function needs to be checked\n\n\n\n\n\n","category":"function"},{"location":"api/#Graphs.savegraph-Tuple{IO, AbstractGraph, SWGFormat}","page":"API reference","title":"Graphs.savegraph","text":"Graphs.savegraph(io, g, SWGFormat())\n\nWrite a graph g with default name \"graph\" to the IO stream io using SWGFormat, and return 1 (number of graphs written).\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.savegraph-Tuple{IO, AbstractGraph, String, SWGFormat}","page":"API reference","title":"Graphs.savegraph","text":"Graphs.savegraph(io, g, gname, SWGFormat())\n\nWrite a graph g with name gname to the IO stream io using SWGFormat, and return 1 (number of graphs written).\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.savegraph-Tuple{IO, Dict, SWGFormat}","page":"API reference","title":"Graphs.savegraph","text":"Graphs.savegraph(io, d, SWGFormat())\n\nWrite a dictionary d of name => graph pairs to the IO stream io using SWGFormat, and return the number of graphs written.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.savegraph-Union{Tuple{U}, Tuple{T}, Tuple{AbstractString, Dict{T, U}}} where {T<:AbstractString, U<:AbstractSimpleWeightedGraph}","page":"API reference","title":"Graphs.savegraph","text":"Graphs.savegraph(fn, d)\n\nwarning: Warning\nThis function needs to be checked\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.weights-Tuple{SimpleWeightedDiGraph}","page":"API reference","title":"Graphs.weights","text":"Graphs.weights(g::SimpleWeightedDiGraph)\n\nReturn the weighted adjacency matrix, stored as an Adjoint.\n\n\n\n\n\n","category":"method"},{"location":"api/#Graphs.weights-Tuple{SimpleWeightedGraph}","page":"API reference","title":"Graphs.weights","text":"Graphs.weights(g::SimpleWeightedGraph)\n\nReturn the weighted adjacency matrix.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs._get_nz_index!-Tuple{SparseArrays.SparseMatrixCSC, Integer, Integer}","page":"API reference","title":"SimpleWeightedGraphs._get_nz_index!","text":"_get_nz_index!(mat::SparseMatrixCSC, i, j)\n\nReturn the index in nzval of mat[i, j]. We assume bounds are already checked\n\nSee https://github.com/JuliaSparse/SparseArrays.jl/blob/fa547689947fadd6c2f3d09ddfcb5f26536f18c8/src/sparsematrix.jl#L2492 for implementation\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.degree_matrix","page":"API reference","title":"SimpleWeightedGraphs.degree_matrix","text":"degree_matrix(g, T; dir)\n\nConstruct the weighted diagonal degree matrix, filled with element type T and considering edge direction dir ∈ [:in, :out, :both] (default is :out).\n\n\n\n\n\n","category":"function"},{"location":"api/#SimpleWeightedGraphs.get_weight-Tuple{AbstractSimpleWeightedGraph, Integer, Integer}","page":"API reference","title":"SimpleWeightedGraphs.get_weight","text":"get_weight(g, u, v)\nget_weight(g, e)\n\nRetrieve the weight of edge (u, v) or e.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.loadswg-Tuple{IO, String}","page":"API reference","title":"SimpleWeightedGraphs.loadswg","text":"loadswg(io, gname)\n\nReturn a single graph with name gname loaded from the IO stream io using SWGFormat.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.loadswg_mult-Tuple{IO}","page":"API reference","title":"SimpleWeightedGraphs.loadswg_mult","text":"loadswg_mult(io)\n\nReturn a dictionary of name => graph pairs loaded from the IO stream io using SWGFormat.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.saveswg-Tuple{IO, AbstractGraph, String}","page":"API reference","title":"SimpleWeightedGraphs.saveswg","text":"saveswg(io, g, gname)\n\nWrite a graph g with name gname to the IO stream io using SWGFormat, and return 1 (number of graphs written).\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.saveswg_mult-Tuple{IO, Dict}","page":"API reference","title":"SimpleWeightedGraphs.saveswg_mult","text":"saveswg_mult(io, d)\n\nWrite a dictionary d of name => graph pairs to the IO stream io using SWGFormat, and return the number of graphs written.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.weight-Tuple{AbstractSimpleWeightedEdge}","page":"API reference","title":"SimpleWeightedGraphs.weight","text":"weight(e)\n\nReturn the weight of a weighted edge.\n\n\n\n\n\n","category":"method"},{"location":"api/#SimpleWeightedGraphs.weighttype-Union{Tuple{AbstractSimpleWeightedGraph{T, U}}, Tuple{T}, Tuple{U}} where {U, T}","page":"API reference","title":"SimpleWeightedGraphs.weighttype","text":"weighttype(g)\n\nReturn the subtype of Real used to represent edge weights.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SimpleWeightedGraphs","category":"page"},{"location":"#SimpleWeightedGraphs","page":"Home","title":"SimpleWeightedGraphs","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SimpleWeightedGraphs.","category":"page"},{"location":"#Quick-start","page":"Home","title":"Quick start","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package defines two new graph types: SimpleWeightedGraph and SimpleWeightedDiGraph. See the tutorial to discover what you can do with them. Also refer to the Graphs.jl package for more complex algorithms.","category":"page"},{"location":"#Caveats","page":"Home","title":"Caveats","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Because SimpleWeighted(Di)Graphs are stored in sparse matrices, they have two major flaws:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Iteratively adding or removing vertices or edges is not very efficient. Building the graph in one go from a list of edge sources, destinations and weights is much faster.\nZero-weight edges are discarded by add_edge!. A possible workaround is to set a very small weight instead.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In additions, self-loops are not supported.","category":"page"},{"location":"#Alternatives","page":"Home","title":"Alternatives","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If your graphs have more than just edge weights to store, take a look at MetaGraphsNext.jl or MetaGraphs.jl for more complex formats.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"CurrentModule = SimpleWeightedGraphs","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> using Graphs, SimpleWeightedGraphs","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Here's how to construct an undirected graph (use SimpleWeightedDiGraph for directed graphs):","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> g = SimpleWeightedGraph(3)\n{3, 0} undirected simple Int64 graph with Float64 weights\n\njulia> add_edge!(g, 1, 2, 0.5);\n\njulia> add_edge!(g, 2, 3, 0.8);\n\njulia> add_edge!(g, 1, 3, 2.0);","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Get the weight of edge from vertex 1 to vertex 2:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> get_weight(g, 1, 2)\n0.5","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Or by passing in any e::AbstractEdge:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> get_weight(g, Edge(1, 2))\n0.5","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Find the shortest path from vertex 1 to vertex 3, taking weights into account:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> enumerate_paths(dijkstra_shortest_paths(g, 1), 3)\n3-element Vector{Int64}:\n 1\n 2\n 3","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Reweight the edge from 1 to 2:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> add_edge!(g, 1, 2, 1.6);","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Rerun the shortest path calculation from 1 to 3:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> enumerate_paths(dijkstra_shortest_paths(g, 1), 3)\n2-element Vector{Int64}:\n 1\n 3","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"It's possible (and faster) to build the graph from arrays of sources, destinations and weights:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> sources = [1, 2, 1];\n\njulia> destinations = [2, 3, 3];\n\njulia> weights = [0.5, 0.8, 2.0];\n\njulia> g = SimpleWeightedGraph(sources, destinations, weights)\n{3, 3} undirected simple Int64 graph with Float64 weights","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The combine keyword handles repeated pairs (sum by default). Unexpected results might occur with non-associative combine functions.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> g = SimpleWeightedGraph([1,2,1], [2,1,2], [1,1,1]; combine = +)\n{2, 1} undirected simple Int64 graph with Int64 weights\n\njulia> g.weights[2,1] == g.weights[1,2] == 3\ntrue","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Notice that weights are indexed by [destination, source] internally:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia> s = SimpleWeightedDiGraph([1,2,1], [2,1,2], [1,1,1]; combine = +);\n\njulia> s.weights[1,2] == 1\ntrue\n\njulia> s.weights[2,1] == 2\ntrue","category":"page"}]
}
