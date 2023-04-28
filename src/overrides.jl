"""
    degree_matrix(g, T; dir)

Construct the weighted diagonal degree matrix, filled with element type `T` and considering edge direction `dir ∈ [:in, :out, :both]` (default is `:out`).
"""
function degree_matrix(
    g::AbstractSimpleWeightedGraph, T::DataType=weighttype(g); dir::Symbol=:out
)
    if is_directed(g)
        if dir == :out
            d = vec(sum(g.weights; dims=1))
        elseif dir == :in
            d = vec(sum(g.weights; dims=2))
        elseif dir == :both
            d = vec(sum(g.weights; dims=1)) + vec(sum(g.weights; dims=2))
        else
            throw(DomainError(dir, "invalid argument, only accept :in, :out and :both"))
        end
    else
        d = vec(sum(g.weights; dims=1))
    end
    return spdiagm(0 => T.(d))
end

"""
    Graphs.adjacency_matrix(g, T; dir)

Construct the weighted adjacency matrix, filled with element type `T` and considering edge direction `dir ∈ [:in, :out, :both]` (default is `:out`).
"""
function Graphs.adjacency_matrix(
    g::AbstractSimpleWeightedGraph, T::DataType=weighttype(g); dir::Symbol=:out
)
    if dir == :out
        return SparseMatrixCSC(T.(copy(g.weights))')
    else
        return T.(copy(g.weights))
    end
end

"""
    Graphs.laplacian_matrix(g, T; dir)

Subtract the adjacency matrix to the degree matrix, both filled with element type `T` and considering edge direction `dir ∈ [:in, :out, :both]` (unlike in Graphs.jl, default is `:out`).
"""
function Graphs.laplacian_matrix(
    g::AbstractSimpleWeightedGraph, T::DataType=weighttype(g); dir::Symbol=:out
)
    return degree_matrix(g, T; dir=dir) - adjacency_matrix(g, T; dir=dir)
end

"""
    Graphs.pagerank(g, α=0.85, n=100, ϵ=1.0e-6)

Apply the page rank algorithm on a weighted graph.
"""
function Graphs.pagerank(g::SimpleWeightedDiGraph, α=0.85, n::Integer=100, ϵ=1.0e-6)
    A = weights(g)
    S = 1 ./ vec(sum(A; dims=2))  # inverse of outdegree
    S[findall(S .== Inf)] .= 0.0
    # scaling the adjmat to stochastic adjacency matrix
    M = (Diagonal(S) * A)'
    N = Int(nv(g))
    # solution vector
    x = fill(1.0 / N, N)
    # personalization vector
    p = fill(1.0 / N, N)
    # temporary to hold the results of SpMV
    y = zeros(Float64, N)
    # adjustment for leaf nodes in digraph
    dangling_weights = p
    is_dangling = findall(S .== 0)
    # save some flops by precomputing this
    pscaled = (1 .- α) .* p
    for _ in 1:n
        xlast = x
        # in place SpMV to conserve memory
        mul!(y, M, x)
        # using broadcast to avoid temporaries
        x = α .* (y .+ sum(x[is_dangling]) .* dangling_weights) .+ pscaled
        # l1 change in solution convergence criterion
        err = sum(abs, (x .- xlast))
        if (err < N * ϵ)
            return x
        end
    end
    return error("Pagerank did not converge after $n iterations.")
end

"""
Graphs.cartesian_product(g, h)

Compute the weighted cartesian product of two weighted graphs.

!!! warning "Warning"
    It is possible that this is suboptimal, but it is the most trivial extension of the implementation used in Graphs.jl.
"""
function Graphs.cartesian_product(g::G, h::G) where {G<:AbstractSimpleWeightedGraph}
    z = G(nv(g) * nv(h))
    id(i, j) = (i - 1) * nv(h) + j
    for e in edges(g)
        i1, i2 = Tuple(e)
        for j in 1:nv(h)
            add_edge!(z, id(i1, j), id(i2, j), weight(e))
        end
    end

    for e in edges(h)
        j1, j2 = Tuple(e)
        for i in vertices(g)
            add_edge!(z, id(i, j1), id(i, j2), weight(e))
        end
    end
    return z
end

# Connected Components on a Sparse Matrix

function _cc(g::SimpleWeightedGraph{T,U}) where {T} where {U}
    a = weights(g)
    comp = 0
    n = size(a, 1)
    marks = zeros(T, n)
    queue = Vector{T}()
    for i in 1:n
        if marks[i] == 0
            comp += 1
            push!(queue, i)
            while !isempty(queue)
                v = pop!(queue)
                marks[v] = comp
                for index in nzrange(a, v)
                    n = a.rowval[index]
                    if marks[n] == 0
                        push!(queue, n)
                    end
                end
            end
        end
    end
    return marks, comp
end

"""
    Graphs.connected_components(g)

Compute the connected components of a weighted graph. Note that an edge with weight `0` will still be counted as an edge if it exists in the sparse weights matrix.
"""
function Graphs.connected_components(g::SimpleWeightedGraph{T,U}) where {T,U}
    marks, num_cc = _cc(g)
    cc = [Vector{T}() for i in 1:num_cc]
    for (i, v) in enumerate(marks)
        push!(cc[v], i)
    end
    return cc
end

"""
    Graphs.induced_subgraph(g, vlist)

Compute the weighted subgraph induced by a list of vertices.

Return a tuple containing the new graph and the list of vertices.
"""
function Graphs.induced_subgraph(
    g::G, vlist::AbstractVector{U}
) where {G<:AbstractSimpleWeightedGraph,U<:Integer}
    T = eltype(g)
    allunique(vlist) || throw(ArgumentError("Vertices in subgraph list must be unique"))
    new_weights = g.weights[T.(vlist), T.(vlist)]
    newg = zero(g)
    newg.weights = new_weights
    return newg, Vector{T}(vlist)
end

function Graphs.induced_subgraph(
    g::G, elist::AbstractVector{E}
) where {G<:AbstractSimpleWeightedGraph} where {E<:AbstractEdge}
    allunique(elist) || throw(ArgumentError("Edges in subgraph list must be unique"))
    T, U = eltype(g), weighttype(g)
    vertex_set = Set{T}()
    for e in elist
        if has_edge(g, e)
            push!(vertex_set, src(e), dst(e))
        else
            @warn "Skipping the edge $(e), since it does not exist in the graph!"
        end
    end
    vertex_list = collect(vertex_set)
    sort!(vertex_list)
    index_map = Dict(vertex_list[i] => i for i in eachindex(vertex_list))
    n = length(vertex_list)
    new_weights = spzeros(weighttype(g), T, n, n)
    I, J, W = T[], T[], U[]
    for e in elist
        if has_edge(g, e)
            i, j = index_map[src(e)], index_map[dst(e)]
            w = get_weight(g, dst(e), src(e))
            push!(I, j)  # storage is transposed!
            push!(J, i)
            push!(W, w)
            if !is_directed(g)
                push!(I, i)
                push!(J, j)
                push!(W, w)
            end
        end
    end
    new_weights = sparse(I, J, W)
    newg = G(new_weights)
    return newg, vertex_list
end
