
"""
    SimpleWeightedGraph{T, U}

A type representing an undirected graph with weights of type `U`.

Note that adding or removing vertices or edges is not particularly performant;
see MetaGraphs.jl for possible alternatives.
"""
mutable struct SimpleWeightedGraph{T<:Integer, U<:Real} <: AbstractSimpleWeightedGraph{T, U}
    weights::SparseMatrixCSC{U,T}
    function SimpleWeightedGraph{T, U}(adjmx::SparseMatrixCSC{U, T}) where {T<:Integer, U<:Real}
        dima,dimb = size(adjmx)
        isequal(dima,dimb) || error("Adjacency / distance matrices must be square")
        issymmetric(adjmx) || error("Adjacency / distance matrices must be symmetric")
        new{T, U}(adjmx)
    end

end

ne(g::SimpleWeightedGraph) = (nnz(g.weights) + nselfloop(g)) ÷ 2

function has_edge(g::SimpleWeightedGraph, u::Integer, v::Integer)
    u ∈ vertices(g) && v ∈ vertices(g) || return false
    _get_nz_index!(g.weights, v, u) != 0 # faster than Base.isstored
end

has_edge(g::SimpleWeightedGraph, e::AbstractEdge) =
    has_edge(g, src(e), dst(e))

function nselfloop(g::SimpleWeightedGraph)
    n = 0
    for i in axes(g.weights, 1)
        n += g.weights[i, i] != 0
    end
    return n
end

SimpleWeightedGraph{T}(adjmx::SparseMatrixCSC{U, T}) where {T <: Integer, U <: Real} =
    SimpleWeightedGraph{T, U}(adjmx)

SimpleWeightedGraph(adjmx::SparseMatrixCSC{U, T}) where {T <: Integer, U <: Real} =
    SimpleWeightedGraph{T, U}(adjmx)

SimpleWeightedGraph(m::AbstractMatrix{U}) where {U <: Real} =
    SimpleWeightedGraph{Int, U}(SparseMatrixCSC{U, Int}(m))
SimpleWeightedGraph{T}(m::AbstractMatrix{U}) where T<:Integer where U<:Real =
    SimpleWeightedGraph{T, U}(SparseMatrixCSC{U, T}(m))
SimpleWeightedGraph{T, U}(m::AbstractMatrix) where T<:Integer where U<:Real =
    SimpleWeightedGraph{T, U}(SparseMatrixCSC{U, T}(m))


SimpleWeightedGraph(g::SimpleWeightedGraph) = SimpleWeightedGraph(copy(g.weights))
function SimpleWeightedGraph{T,U}(g::SimpleWeightedGraph) where {T<:Integer, U<:Real}
    return SimpleWeightedGraph(SparseMatrixCSC{U, T}(copy(g.weights)))
end

# Graph{UInt8}(6), Graph{Int16}(7), Graph{UInt8}()
function (::Type{SimpleWeightedGraph{T, U}})(n::Integer = 0) where T<:Integer where U<:Real
    weights = spzeros(U, T, T(n), T(n))
    return SimpleWeightedGraph{T, U}(weights)
end

# Graph()
SimpleWeightedGraph() = SimpleWeightedGraph(Matrix{Float64}(undef, 0, 0))

# Graph(6), Graph(0x5)
SimpleWeightedGraph(n::T) where T<:Integer = SimpleWeightedGraph{T, Float64}(n)

# Graph(UInt8)
SimpleWeightedGraph(::Type{T}) where T<:Integer = SimpleWeightedGraph{T, Float64}(zero(T))

# Graph(UInt8, Float32)
SimpleWeightedGraph(::Type{T}, ::Type{U}) where {T<:Integer, U<:Real} = SimpleWeightedGraph{T, U}(zero(T))

# Graph(SimpleGraph)

function SimpleWeightedGraph(g::Graphs.AbstractGraph{T}, ::Type{U}=Float64) where {T <: Integer, U <: Real}
    adj_matrix = if Graphs.is_directed(g)
        # TODO abstract function instead of SimpleGraph constructor
        adjacency_matrix(Graphs.SimpleGraphs.SimpleGraph(g), U)
    else
        adjacency_matrix(g, U)
    end
    return SimpleWeightedGraph{T, U}(adj_matrix)
end

function SimpleWeightedGraph(g::Graphs.AbstractGraph{T}, x::U) where {T <: Integer, U <: Real}
    adj_matrix = if Graphs.is_directed(g)
        # TODO abstract function instead of SimpleGraph constructor
        adjacency_matrix(Graphs.SimpleGraphs.SimpleGraph(g), U)
    else
        adjacency_matrix(g, U)
    end
    return SimpleWeightedGraph{T, U}(x .* adj_matrix)
end

# SimpleWeightedGraph{T, U}(SimpleGraph)
function (::Type{SimpleWeightedGraph{T, U}})(g::Graphs.AbstractGraph)  where {T<:Integer, U <: Real}
    adj_matrix = if Graphs.is_directed(g)
        # TODO abstract function instead of SimpleGraph constructor
        adjacency_matrix(Graphs.SimpleGraphs.SimpleGraph{T}(g), U)
    else
        adjacency_matrix(g, U)
    end
    return SimpleWeightedGraph{T, U}(adj_matrix)
end

# Graph(srcs, dsts, weights)
function SimpleWeightedGraph(i::AbstractVector{T}, j::AbstractVector{T}, v::AbstractVector{U}; combine = +) where {T<:Integer, U<:Real}
    m = max(maximum(i), maximum(j))
    s = sparse(vcat(i,j), vcat(j,i), vcat(v,v), m, m, combine)
    SimpleWeightedGraph{T, U}(s)
end

Graphs.SimpleGraph(g::SimpleWeightedGraph) = SimpleGraph(g.weights)

edgetype(::SimpleWeightedGraph{T, U}) where {T<:Integer, U<:Real} = SimpleWeightedGraphEdge{T,U}

edges(g::SimpleWeightedGraph) = (SimpleWeightedEdge(x[1], x[2], x[3]) for x in zip(findnz(triu(g.weights))...))
weights(g::SimpleWeightedGraph) = g.weights

function outneighbors(g::SimpleWeightedGraph, v::Integer)
    mat = g.weights
    return view(mat.rowval, mat.colptr[v]:mat.colptr[v+1]-1)
end
inneighbors(g::SimpleWeightedGraph, x...) = outneighbors(g, x...)

# add_edge! will overwrite weights.
function add_edge!(g::SimpleWeightedGraph, e::SimpleWeightedGraphEdge)
    T = eltype(g)
    U = weighttype(g)
    s_, d_, w = Tuple(e)

    if w == zero(U)
        @warn "Note: adding edges with a zero weight to this graph type has no effect." maxlog=1 _id=:swg_add_edge_zero
        return false
    end

    s = T(s_)
    d = T(d_)
    (s in vertices(g) && d in vertices(g)) || return false
    @inbounds g.weights[d, s] = w
    @inbounds g.weights[s, d] = w
    return true
end

rem_edge!(g::SimpleWeightedGraph{T, U}, e::AbstractEdge) where {T<:Integer, U<:Real} =
    rem_edge!(g, src(e), dst(e))

function rem_edge!(g::SimpleWeightedGraph{T, U}, u::Integer, v::Integer) where {T<:Integer, U<:Real}
    u ∈ vertices(g) && v ∈ vertices(g) || return false
    w = g.weights
    indx = _get_nz_index!(w, u, v)
    indx == 0 && return false
    @view(w.colptr[v+one(v):end]) .-= T(1)
    deleteat!(w.rowval, indx)
    deleteat!(w.nzval, indx)
    (u == v) && return true
    indx = _get_nz_index!(w, v, u)
    @view(w.colptr[u+one(u):end]) .-= T(1)
    deleteat!(w.rowval, indx)
    deleteat!(w.nzval, indx)
    return true
end

@doc_str """
    rem_vertex!(g::SimpleWeightedGraph, v)

Remove the vertex `v` from graph `g`. Return false if removal fails
(e.g., if vertex is not in the graph); true otherwise.

### Implementation Notes
This operation has to be performed carefully if one keeps external
data structures indexed by edges or vertices in the graph, since
internally the removal results in all vertices with indices greater than `v`
being shifted down one.
"""
function rem_vertex!(g::SimpleWeightedGraph, v::Integer)
    v in vertices(g) || return false
    n = nv(g)
    g.weights = g.weights[1:n .!= v, 1:n .!= v]
    return true
end


==(g::SimpleWeightedGraph, h::SimpleWeightedGraph) = g.weights == h.weights

is_directed(::Type{<:SimpleWeightedGraph}) = false

function Base.getindex(g::SimpleWeightedGraph{T, U}, e::AbstractEdge, ::Val{:weight}) where {T, U, S}
    return g.weights[src(e), dst(e)]
end

function Base.getindex(g::SimpleWeightedGraph{T, U}, i::Integer, j::Integer, ::Val{:weight}) where {T, U, S}
    return g.weights[i, j]
end
