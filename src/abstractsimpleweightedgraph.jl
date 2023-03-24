
"""
    AbstractSimpleWeightedGraph

An abstract type representing a simple graph structure.

# Required fields
- `weightmx::AbstractSparseMatrix{Real}`
"""
abstract type AbstractSimpleWeightedGraph{T<:Integer,U<:Real} <: AbstractGraph{T} end

function Base.show(io::IO, g::AbstractSimpleWeightedGraph{T,U}) where {T} where {U}
    dir = is_directed(g) ? "directed" : "undirected"
    return print(io, "{$(nv(g)), $(ne(g))} $dir simple $T graph with $U weights")
end

## Interface

Base.eltype(::AbstractSimpleWeightedGraph{T,U}) where {T} where {U} = T
Base.zero(::T) where {T<:AbstractSimpleWeightedGraph} = T()
Base.copy(g::T) where {T<:AbstractSimpleWeightedGraph} = T(copy(weights(g)))

weighttype(::AbstractSimpleWeightedGraph{T,U}) where {T} where {U} = U
get_weight(g::AbstractSimpleWeightedGraph, u::Integer, v::Integer) = weights(g)[v, u]

## Vertices

Graphs.nv(g::AbstractSimpleWeightedGraph{T,U}) where {T} where {U} = T(size(weights(g), 1))
Graphs.vertices(g::AbstractSimpleWeightedGraph{T,U}) where {T} where {U} = one(T):nv(g)
Graphs.has_vertex(g::AbstractSimpleWeightedGraph, v::Integer) = v in vertices(g)

# TODO: manipulate SparseMatrixCSC directly
Graphs.add_vertex!(g::AbstractSimpleWeightedGraph) = add_vertices!(g, 1)

## Edges

# handles single-argument edge constructors such as pairs and tuples
function Graphs.has_edge(g::AbstractSimpleWeightedGraph{T,U}, x) where {T} where {U}
    return has_edge(g, edgetype(g)(x))
end

function Graphs.add_edge!(g::AbstractSimpleWeightedGraph{T,U}, x) where {T} where {U}
    return add_edge!(g, edgetype(g)(x))
end

# handles two-argument edge constructors like src,dst
Graphs.has_edge(g::AbstractSimpleWeightedGraph, x, y) = has_edge(g, edgetype(g)(x, y, 0))
Graphs.add_edge!(g::AbstractSimpleWeightedGraph, x, y) = add_edge!(g, edgetype(g)(x, y, 1))

function Graphs.add_edge!(g::AbstractSimpleWeightedGraph, x, y, z)
    return add_edge!(g, edgetype(g)(x, y, z))
end

function Graphs.rem_edge!(
    g::AbstractSimpleWeightedGraph{T,U}, u::Integer, v::Integer
) where {T,U}
    return rem_edge!(g, edgetype(g)(T(u), T(v), one(U)))
end

## Miscellaneous

function Base.issubset(g::T, h::T) where {T<:AbstractSimpleWeightedGraph}
    (gmin, gmax) = extrema(vertices(g))
    (hmin, hmax) = extrema(vertices(h))
    return (hmin <= gmin <= gmax <= hmax) && issubset(edges(g), edges(h))
end
