module SimpleWeightedGraphs

using LinearAlgebra
using Markdown
using SparseArrays: SparseMatrixCSC, sparse, spzeros, nnz, findnz, spdiagm, nzrange

using Graphs: Graphs
using Graphs: AbstractGraph, AbstractEdge, AbstractEdgeIter, AbstractGraphFormat
using Graphs: SimpleGraph, SimpleDiGraph
using Graphs: src, dst
using Graphs: edgetype, is_directed, nv, ne, vertices, edges
using Graphs: add_vertex!, add_vertices!, add_edge!, rem_vertex!, rem_edge!
using Graphs: has_vertex, has_edge, inneighbors, outneighbors
using Graphs: indegree, outdegree, degree, has_self_loops, num_self_loops
using Graphs: adjacency_matrix, laplacian_matrix, weights
using Graphs: connected_components, cartesian_product, induced_subgraph, pagerank
using Graphs: loadgraph, loadgraphs, savegraph
using Graphs: _NI

export AbstractSimpleWeightedGraph, AbstractSimpleWeightedEdge
export SimpleWeightedGraph, SimpleWeightedDiGraph
export SimpleWeightedEdge, SimpleWeightedGraphEdge, SimpleWeightedDiGraphEdge
export WGraph, WDiGraph, SWGFormat
export weight, weighttype, get_weight, degree_matrix

include("utils.jl")
include("simpleweightededge.jl")
include("abstractsimpleweightedgraph.jl")
include("simpleweighteddigraph.jl")
include("simpleweightedgraph.jl")
include("conversion.jl")
include("overrides.jl")
include("persistence.jl")

end
