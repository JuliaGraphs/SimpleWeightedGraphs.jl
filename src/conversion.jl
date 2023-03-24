SimpleWeightedDiGraph(g::SimpleWeightedGraph) = SimpleWeightedDiGraph(copy(g.weights))

function SimpleWeightedDiGraph{T,U}(g::SimpleWeightedGraph) where {T<:Integer,U<:Real}
    return SimpleWeightedDiGraph(SparseMatrixCSC{U,T}(copy(g.weights)))
end

SimpleWeightedGraph(g::SimpleWeightedDiGraph) = SimpleWeightedGraph(g.weights .+ g.weights')

function SimpleWeightedGraph{T,U}(g::SimpleWeightedDiGraph) where {T<:Integer,U<:Real}
    return SimpleWeightedGraph(SparseMatrixCSC{U,T}(g.weights .+ g.weights'))
end

function Base.convert(
    ::Type{SparseMatrixCSC{U,T}}, g::AbstractSimpleWeightedGraph
) where {U<:Real,T<:Integer}
    return SparseMatrixCSC{U,T}(g.weights)
end
