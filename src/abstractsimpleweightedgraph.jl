
"""
    AbstractSimpleWeightedGraph

An abstract type representing a simple graph structure with edge weights.

The only requirement for concrete subtypes is that they should implement a method `Graphs.weights(g)` which returns the weighted adjacency matrix.
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

"""
    weighttype(g)

Return the subtype of `Real` used to represent edge weights.
"""
weighttype(::AbstractSimpleWeightedGraph{T,U}) where {T} where {U} = U

"""
    get_weight(g, u, v)

Retrieve the weight of edge `(u, v)`.
"""
get_weight(g::AbstractSimpleWeightedGraph, u::Integer, v::Integer) = weights(g)[u, v]

## Vertices

Graphs.nv(g::AbstractSimpleWeightedGraph{T,U}) where {T} where {U} = T(size(weights(g), 1))
Graphs.vertices(g::AbstractSimpleWeightedGraph{T,U}) where {T} where {U} = one(T):nv(g)
Graphs.has_vertex(g::AbstractSimpleWeightedGraph, v::Integer) = v in vertices(g)

# TODO: manipulate SparseMatrixCSC directly
Graphs.add_vertex!(g::AbstractSimpleWeightedGraph) = add_vertices!(g, 1)

function Graphs.add_vertices!(g::AbstractSimpleWeightedGraph, n::Integer)
    T = eltype(g)
    U = weighttype(g)
    (nv(g) + one(T) <= nv(g)) && return false       # test for overflow
    emptycols = spzeros(U, nv(g) + n, n)
    g.weights = hcat(g.weights, emptycols[1:nv(g), :])
    g.weights = vcat(g.weights, emptycols')
    return true
end

## Edges

# Handle single-argument edge constructors such as pairs and tuples

"""
    Graphs.has_edge(g, e)

Check the existence of the edge `e` in the graph, where `e` is any object that can be converted into `edgetype(g)`.
"""
function Graphs.has_edge(g::AbstractSimpleWeightedGraph{T,U}, e) where {T,U}
    return has_edge(g, edgetype(g)(e))
end

"""
    Graphs.add_edge!(g, e)

Add the edge `e` to the graph, where `e` is any object that can be converted into `edgetype(g)`.
"""
function Graphs.add_edge!(g::AbstractSimpleWeightedGraph{T,U}, e) where {T,U}
    return add_edge!(g, edgetype(g)(e))
end

# Handle two-argument edge constructors like src,dst

"""
    Graphs.has_edge(g, u, v)

Check the existence of the edge `(u, v)` in the graph.
"""
function Graphs.has_edge(g::AbstractSimpleWeightedGraph{T,U}, u, v) where {T,U}
    return has_edge(g, edgetype(g)(u, v, 0))
end

"""
    Graphs.add_edge!(g, u, v)

Add the edge `(u, v)` to the graph with a default weight of 1.
"""
function Graphs.add_edge!(g::AbstractSimpleWeightedGraph{T,U}, u, v) where {T,U}
    return add_edge!(g, edgetype(g)(T(u), T(v), one(U)))
end

"""
    Graphs.rem_edge!(g, u, v)

Remove the edge `(u, v)` from the graph.
"""
function Graphs.rem_edge!(
    g::AbstractSimpleWeightedGraph{T,U}, u::Integer, v::Integer
) where {T,U}
    return rem_edge!(g, edgetype(g)(T(u), T(v), one(U)))
end

# Handle three-argument constructors

"""
    Graphs.add_edge!(g, u, v, w)

Add the edge `(u, v)` to the graph with a weight of `w`.
"""
function Graphs.add_edge!(g::AbstractSimpleWeightedGraph{T,U}, u, v, w) where {T,U}
    return add_edge!(g, edgetype(g)(T(u), T(v), U(w)))
end

## Miscellaneous

function Base.issubset(g::T, h::T) where {T<:AbstractSimpleWeightedGraph}
    (gmin, gmax) = extrema(vertices(g))
    (hmin, hmax) = extrema(vertices(h))
    return (hmin <= gmin <= gmax <= hmax) && issubset(edges(g), edges(h))
end
