
"""
    SimpleWeightedGraph{T, U}

A type representing an undirected weighted graph with vertices of type `T` and edge weights of type `U`.

# Fields
- `weights::SparseMatrixCSC{U,T}`: weighted adjacency matrix, indexed by `(dst, src)`

!!! tip "Performance"
    Iteratively adding/removing vertices or edges is not very efficient for this type of graph: better construct the graph in one shot if possible.

# Basic constructors
```
SimpleWeightedGraph()  # empty
SimpleWeightedGraph(n)  # n vertices, no edges
SimpleWeightedGraph(graph)  # from graph
SimpleWeightedGraph(adjmx)  # from adjacency matrix
SimpleWeightedGraph(sources, destinations, weights)  # from list of edges
```
Use `methods(SimpleWeightedGraph)` for the full list of constructors.
"""
mutable struct SimpleWeightedGraph{T<:Integer,U<:Real} <: AbstractSimpleWeightedGraph{T,U}
    weights::SparseMatrixCSC{U,T}
    function SimpleWeightedGraph{T,U}(
        adjmx::SparseMatrixCSC{U,T}
    ) where {T<:Integer,U<:Real}
        dima, dimb = size(adjmx)
        isequal(dima, dimb) || error("Adjacency / distance matrices must be square")
        issymmetric(adjmx) || error("Adjacency / distance matrices must be symmetric")
        return new{T,U}(adjmx)
    end
end

"""
    WGraph

Alias for `SimpleWeightedGraph`.
"""
const WGraph = SimpleWeightedGraph

Graphs.ne(g::SimpleWeightedGraph) = (nnz(g.weights) + nselfloop(g)) ÷ 2

function Graphs.has_edge(g::SimpleWeightedGraph, u::Integer, v::Integer)
    u in vertices(g) && v in vertices(g) || return false
    return _get_nz_index!(g.weights, v, u) != 0 # faster than Base.isstored
end

Graphs.has_edge(g::SimpleWeightedGraph, e::AbstractEdge) = has_edge(g, src(e), dst(e))

function nselfloop(g::SimpleWeightedGraph)
    n = 0
    for i in axes(g.weights, 1)
        n += g.weights[i, i] != 0
    end
    return n
end

function SimpleWeightedGraph{T}(adjmx::SparseMatrixCSC{U,T}) where {T<:Integer,U<:Real}
    return SimpleWeightedGraph{T,U}(adjmx)
end

function SimpleWeightedGraph(adjmx::SparseMatrixCSC{U,T}) where {T<:Integer,U<:Real}
    return SimpleWeightedGraph{T,U}(adjmx)
end

function SimpleWeightedGraph(m::AbstractMatrix{U}) where {U<:Real}
    return SimpleWeightedGraph{Int,U}(SparseMatrixCSC{U,Int}(m))
end

function SimpleWeightedGraph{T}(m::AbstractMatrix{U}) where {T<:Integer,U<:Real}
    return SimpleWeightedGraph{T,U}(SparseMatrixCSC{U,T}(m))
end

function SimpleWeightedGraph{T,U}(m::AbstractMatrix) where {T<:Integer,U<:Real}
    return SimpleWeightedGraph{T,U}(SparseMatrixCSC{U,T}(m))
end

SimpleWeightedGraph(g::SimpleWeightedGraph) = SimpleWeightedGraph(copy(g.weights))
function SimpleWeightedGraph{T,U}(g::SimpleWeightedGraph) where {T<:Integer,U<:Real}
    return SimpleWeightedGraph(SparseMatrixCSC{U,T}(copy(g.weights)))
end

# Graph{UInt8}(6), Graph{Int16}(7), Graph{UInt8}()
function (::Type{SimpleWeightedGraph{T,U}})(n::Integer=0) where {T<:Integer} where {U<:Real}
    weights = spzeros(U, T, T(n), T(n))
    return SimpleWeightedGraph{T,U}(weights)
end

# Graph()
SimpleWeightedGraph() = SimpleWeightedGraph(Matrix{Float64}(undef, 0, 0))

# Graph(6), Graph(0x5)
SimpleWeightedGraph(n::T) where {T<:Integer} = SimpleWeightedGraph{T,Float64}(n)

# Graph(UInt8)
SimpleWeightedGraph(::Type{T}) where {T<:Integer} = SimpleWeightedGraph{T,Float64}(zero(T))

# Graph(UInt8, Float32)
function SimpleWeightedGraph(::Type{T}, ::Type{U}) where {T<:Integer,U<:Real}
    return SimpleWeightedGraph{T,U}(zero(T))
end

# Graph(SimpleGraph)

function SimpleWeightedGraph(
    g::Graphs.AbstractGraph{T}, ::Type{U}=Float64
) where {T<:Integer,U<:Real}
    adj_matrix = if is_directed(g)
        # TODO: abstract function instead of SimpleGraph constructor
        adjacency_matrix(SimpleGraph(g), U)
    else
        adjacency_matrix(g, U)
    end
    return SimpleWeightedGraph{T,U}(adj_matrix)
end

function SimpleWeightedGraph(g::Graphs.AbstractGraph{T}, x::U) where {T<:Integer,U<:Real}
    adj_matrix = if is_directed(g)
        # TODO: abstract function instead of SimpleGraph constructor
        adjacency_matrix(SimpleGraph(g), U)
    else
        adjacency_matrix(g, U)
    end
    return SimpleWeightedGraph{T,U}(x .* adj_matrix)
end

# SimpleWeightedGraph{T, U}(SimpleGraph)
function (::Type{SimpleWeightedGraph{T,U}})(
    g::Graphs.AbstractGraph
) where {T<:Integer,U<:Real}
    adj_matrix = if is_directed(g)
        # TODO abstract function instead of SimpleGraph constructor
        adjacency_matrix(SimpleGraph{T}(g), U)
    else
        adjacency_matrix(g, U)
    end
    return SimpleWeightedGraph{T,U}(adj_matrix)
end

# Graph(srcs, dsts, weights)
function SimpleWeightedGraph(
    i::AbstractVector{T}, j::AbstractVector{T}, v::AbstractVector{U}; combine=+
) where {T<:Integer,U<:Real}
    m = max(maximum(i), maximum(j))
    s = sparse(vcat(i, j), vcat(j, i), vcat(v, v), m, m, combine)
    return SimpleWeightedGraph{T,U}(s)
end

Graphs.SimpleGraph(g::SimpleWeightedGraph) = SimpleGraph(g.weights)

function Graphs.edgetype(::SimpleWeightedGraph{T,U}) where {T<:Integer,U<:Real}
    return SimpleWeightedGraphEdge{T,U}
end

function Graphs.edges(g::SimpleWeightedGraph)
    return (SimpleWeightedEdge(x[1], x[2], x[3]) for x in zip(findnz(triu(g.weights))...))
end

"""
    Graphs.weights(g::SimpleWeightedGraph)

Return the weighted adjacency matrix.
"""
Graphs.weights(g::SimpleWeightedGraph) = g.weights

function Graphs.outneighbors(g::SimpleWeightedGraph, v::Integer)
    mat = g.weights
    return view(mat.rowval, mat.colptr[v]:(mat.colptr[v + 1] - 1))
end

Graphs.inneighbors(g::SimpleWeightedGraph, x...) = outneighbors(g, x...)

# add_edge! will overwrite weights.
function Graphs.add_edge!(g::SimpleWeightedGraph, e::SimpleWeightedGraphEdge)
    T = eltype(g)
    U = weighttype(g)
    s_, d_, w = Tuple(e)

    if w == zero(U)
        @warn "Note: adding edges with a zero weight to this graph type has no effect." maxlog =
            1 _id = :swg_add_edge_zero
        return false
    end

    s = T(s_)
    d = T(d_)
    (s in vertices(g) && d in vertices(g)) || return false
    @inbounds g.weights[d, s] = w
    @inbounds g.weights[s, d] = w
    return true
end

"""
    Graphs.rem_edge!(g, e)

Remove the edge `e` from the graph.
"""
Graphs.rem_edge!(g::SimpleWeightedGraph, e::AbstractEdge) = rem_edge!(g, src(e), dst(e))

function Graphs.rem_edge!(
    g::SimpleWeightedGraph{T,U}, u::Integer, v::Integer
) where {T<:Integer,U<:Real}
    (u in vertices(g) && v in vertices(g)) || return false
    w = g.weights
    indx_uv = _get_nz_index!(w, u, v) # get the index in nzval
    indx_uv == 0 && return false # the edge does not exist
    @view(w.colptr[(v + one(v)):end]) .-= T(1) # there is one value less in column v
    # we remove the stored value
    deleteat!(w.rowval, indx_uv)
    deleteat!(w.nzval, indx_uv)
    (u == v) && return true
    # same for the reverse edge
    indx_vu = _get_nz_index!(w, v, u)
    @view(w.colptr[(u + one(u)):end]) .-= T(1)
    deleteat!(w.rowval, indx_vu)
    deleteat!(w.nzval, indx_vu)
    return true
end

"""
    rem_vertex!(g::SimpleWeightedGraph, v)

Remove the vertex `v` from graph `g`. Return false if removal fails (e.g., if vertex is not in the graph) and true otherwise.

!!! tip "Correctness"
    This operation has to be performed carefully if one keeps external data structures indexed by edges or vertices in the graph, since internally the removal results in all vertices with indices greater than `v` being shifted down one.
"""
function Graphs.rem_vertex!(g::SimpleWeightedGraph, v::Integer)
    v in vertices(g) || return false
    n = nv(g)
    all_except_v = (1:n) .!= v
    g.weights = g.weights[all_except_v, all_except_v]
    return true
end

Base.:(==)(g::SimpleWeightedGraph, h::SimpleWeightedGraph) = g.weights == h.weights

Graphs.is_directed(::Type{<:SimpleWeightedGraph}) = false

"""
    g[e, :weight]

Return the weight of edge `e`.
"""
function Base.getindex(g::SimpleWeightedGraph, e::AbstractEdge, ::Val{:weight})
    return g.weights[src(e), dst(e)]
end

"""
    g[i, j, :weight]

Return the weight of edge `(i, j)`.
"""
function Base.getindex(g::SimpleWeightedGraph, i::Integer, j::Integer, ::Val{:weight})
    return g.weights[i, j]
end
