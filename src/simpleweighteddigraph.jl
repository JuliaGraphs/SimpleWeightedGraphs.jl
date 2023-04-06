"""
    SimpleWeightedDiGraph{T,U}

A type representing a directed weighted graph with vertices of type `T` and edge weights of type `U`.

# Fields
- `weights::SparseMatrixCSC{U,T}`: weighted adjacency matrix, indexed by `(dst, src)`

!!! tip "Performance"
    Iteratively adding/removing vertices or edges is not very efficient for this type of graph: better construct the graph in one shot if possible.

# Basic constructors
```
SimpleWeightedDiGraph()  # empty
SimpleWeightedDiGraph(n)  # n vertices, no edges
SimpleWeightedDiGraph(graph)  # from graph
SimpleWeightedDiGraph(adjmx; permute)  # from adjacency matrix, possibly transposed
SimpleWeightedDiGraph(sources, destinations, weights)  # from list of edges
```
Use `methods(SimpleWeightedDiGraph)` for the full list of constructors.
"""
mutable struct SimpleWeightedDiGraph{T<:Integer,U<:Real} <: AbstractSimpleWeightedGraph{T,U}
    weights::SparseMatrixCSC{U,T}
    function SimpleWeightedDiGraph{T,U}(
        adjmx::SparseMatrixCSC{U,T}; permute=true
    ) where {T<:Integer,U<:Real}
        dima, dimb = size(adjmx)
        isequal(dima, dimb) || error("Adjacency / distance matrices must be square")
        return permute ? new{T,U}(permutedims(adjmx)) : new{T,U}(adjmx)
    end
end

"""
    WDiGraph

Alias for `SimpleWeightedDiGraph`.
"""
const WDiGraph = SimpleWeightedDiGraph

function SimpleWeightedDiGraph{T}(
    adjmx::SparseMatrixCSC{U,T}; permute=true
) where {T<:Integer,U<:Real}
    return SimpleWeightedDiGraph{T,U}(adjmx; permute=permute)
end

function SimpleWeightedDiGraph(
    adjmx::SparseMatrixCSC{U,T}; permute=true
) where {T<:Integer,U<:Real}
    return SimpleWeightedDiGraph{T,U}(adjmx; permute=permute)
end

function SimpleWeightedDiGraph(m::AbstractMatrix{U}) where {U<:Real}
    return SimpleWeightedDiGraph{Int,U}(SparseMatrixCSC{U,Int}(m))
end

function SimpleWeightedDiGraph{T}(m::AbstractMatrix{U}) where {T<:Integer,U<:Real}
    return SimpleWeightedDiGraph{T,U}(SparseMatrixCSC{U,T}(m))
end

function SimpleWeightedDiGraph{T,U}(m::AbstractMatrix) where {T<:Integer,U<:Real}
    return SimpleWeightedDiGraph{T,U}(SparseMatrixCSC{U,T}(m))
end

function SimpleWeightedDiGraph(g::SimpleWeightedDiGraph)
    return SimpleWeightedDiGraph(copy(g.weights); permute=false)
end

function SimpleWeightedDiGraph{T,U}(g::SimpleWeightedDiGraph) where {T<:Integer,U<:Real}
    return SimpleWeightedDiGraph(SparseMatrixCSC{U,T}(copy(g.weights)); permute=false)
end

Graphs.ne(g::SimpleWeightedDiGraph) = nnz(g.weights)

function Graphs.has_edge(g::SimpleWeightedDiGraph, u::Integer, v::Integer)
    (u in vertices(g) && v in vertices(g)) || return false
    return _get_nz_index!(g.weights, v, u) != 0 # faster than Base.isstored
end

Graphs.has_edge(g::SimpleWeightedDiGraph, e::AbstractEdge) = has_edge(g, src(e), dst(e))

function SimpleWeightedDiGraph{T,U}(n::Integer=0) where {T<:Integer,U<:Real}
    weights = spzeros(U, T, T(n), T(n))
    return SimpleWeightedDiGraph{T,U}(weights)
end

# Graph()
SimpleWeightedDiGraph() = SimpleWeightedDiGraph{Int,Float64}()

# Graph(6), Graph(0x5)
SimpleWeightedDiGraph(n::T) where {T<:Integer} = SimpleWeightedDiGraph{T,Float64}(n)

# Graph(UInt8)
function SimpleWeightedDiGraph(::Type{T}) where {T<:Integer}
    return SimpleWeightedDiGraph{T,Float64}(zero(T))
end

# Graph(UInt8, Float32)
function SimpleWeightedDiGraph(::Type{T}, ::Type{U}) where {T<:Integer,U<:Real}
    return SimpleWeightedDiGraph{T,U}(zero(T))
end

# DiGraph(AbstractGraph, ::Type{U})
function SimpleWeightedDiGraph(
    g::Graphs.AbstractGraph{T}, ::Type{U}=Float64
) where {U<:Real,T}
    return SimpleWeightedDiGraph{T}(adjacency_matrix(g, U))
end

function SimpleWeightedDiGraph(g::Graphs.AbstractGraph{T}, x::U) where {U<:Real,T}
    m = adjacency_matrix(g, U)'
    return SimpleWeightedDiGraph{T,U}(x .* m; permute=false)
end

# DiGraph(srcs, dsts, weights)
function SimpleWeightedDiGraph(
    i::AbstractVector{T}, j::AbstractVector{T}, v::AbstractVector{U}; combine=+
) where {T<:Integer,U<:Real}
    m = max(maximum(j), maximum(i))
    return SimpleWeightedDiGraph{T,U}(sparse(j, i, v, m, m, combine); permute=false)
end

Graphs.SimpleDiGraph(g::SimpleWeightedDiGraph) = SimpleDiGraph(g.weights')

function Graphs.edgetype(::SimpleWeightedDiGraph{T,U}) where {T<:Integer} where {U<:Real}
    return SimpleWeightedGraphEdge{T,U}
end

function Graphs.edges(g::SimpleWeightedDiGraph)
    return (SimpleWeightedEdge(x[2], x[1], x[3]) for x in zip(findnz(g.weights)...))
end

"""
    Graphs.weights(g::SimpleWeightedDiGraph)

Return the weighted adjacency matrix, stored as an `Adjoint`.
"""
Graphs.weights(g::SimpleWeightedDiGraph) = g.weights'

"""
    Graphs.outneighbors(g::SimpleWeightedDiGraph, v)

Return the vector of outneighbors of vertex `v`.

!!! tip "Performance"
    This function is more efficient than `inneighbors` for directed weighted graphs.
"""
function Graphs.outneighbors(g::SimpleWeightedDiGraph, v::Integer)
    mat = g.weights
    return view(mat.rowval, mat.colptr[v]:(mat.colptr[v + 1] - 1))
end

"""
    Graphs.inneighbors(g::SimpleWeightedDiGraph, v)

Return the vector of inneighbors of vertex `v`.

!!! tip "Performance"
    This function is less efficient than `inneighbors` for directed weighted graphs (it allocates a new vector).
"""
Graphs.inneighbors(g::SimpleWeightedDiGraph, v::Integer) = g.weights[v, :].nzind

# add_edge! will overwrite weights.
function Graphs.add_edge!(g::SimpleWeightedDiGraph, e::SimpleWeightedGraphEdge)
    T = eltype(g)
    U = weighttype(g)
    s_, d_, w = Tuple(e)

    if w == zero(U)
        @warn "Note: adding edges with a zero weight to this graph type has no effect." maxlog =
            1 _id = :swd_add_edge_zero
        return false
    end

    s = T(s_)
    d = T(d_)
    (s in vertices(g) && d in vertices(g)) || return false
    @inbounds g.weights[d, s] = w
    return true
end

"""
    Graphs.rem_edge!(g, e)

Remove the edge `e` from the graph.
"""
Graphs.rem_edge!(g::SimpleWeightedDiGraph, e::AbstractEdge) = rem_edge!(g, src(e), dst(e))

function Graphs.rem_edge!(g::SimpleWeightedDiGraph{T}, u::Integer, v::Integer) where {T}
    (u in vertices(g) && v in vertices(g)) || return false
    w = g.weights
    indx = _get_nz_index!(w, v, u) # get the index in nzval
    indx == 0 && return false # the edge does not exist
    @view(w.colptr[(u + one(u)):end]) .-= T(1) # there is one value less in column u
    # we remove the stored value
    deleteat!(w.rowval, indx)
    deleteat!(w.nzval, indx)

    return true
end

"""
    rem_vertex!(g::SimpleWeightedDiGraph, v)

Remove the vertex `v` from graph `g`. Return false if removal fails (e.g., if vertex is not in the graph) and true otherwise.

!!! tip "Correctness"
    This operation has to be performed carefully if one keeps external data structures indexed by edges or vertices in the graph, since internally the removal results in all vertices with indices greater than `v` being shifted down one.
"""
function Graphs.rem_vertex!(g::SimpleWeightedDiGraph, v::Integer)
    v in vertices(g) || return false
    n = nv(g)
    all_except_v = (1:n) .!= v
    g.weights = g.weights[all_except_v, all_except_v]
    return true
end

Base.copy(g::SimpleWeightedDiGraph) = SimpleWeightedDiGraph(copy(g.weights'))

Base.:(==)(g::SimpleWeightedDiGraph, h::SimpleWeightedDiGraph) = g.weights == h.weights

Graphs.is_directed(::Type{<:SimpleWeightedDiGraph}) = true

"""
    g[e, Val(:weight)]

Return the weight of edge `e`.
"""
function Base.getindex(g::SimpleWeightedDiGraph, e::AbstractEdge, ::Val{:weight})
    return g.weights[dst(e), src(e)]
end

"""
    g[i, j, Val(:weight)]

Return the weight of edge `(i, j)`.
"""
function Base.getindex(g::SimpleWeightedDiGraph, i::Integer, j::Integer, ::Val{:weight})
    return g.weights[j, i]
end
