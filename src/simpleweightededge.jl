abstract type AbstractSimpleWeightedEdge{T} <: AbstractEdge{T} end

struct SimpleWeightedEdge{T<:Integer,U<:Real} <: AbstractSimpleWeightedEdge{T}
    src::T
    dst::T
    weight::U
end

const SimpleWeightedGraphEdge = SimpleWeightedEdge
const SimpleWeightedDiGraphEdge = SimpleWeightedEdge

SimpleWeightedEdge(t::NTuple{2}) = SimpleWeightedEdge(t[1], t[2], one(Float64))
SimpleWeightedEdge(t::NTuple{3}) = SimpleWeightedEdge(t[1], t[2], t[3])
SimpleWeightedEdge(p::Pair) = SimpleWeightedEdge(p.first, p.second, one(Float64))

function SimpleWeightedEdge{T,U}(p::Pair) where {T<:Integer} where {U<:Real}
    return SimpleWeightedEdge(T(p.first), T(p.second), one(U))
end

function SimpleWeightedEdge{T,U}(t::NTuple{3}) where {T<:Integer} where {U<:Real}
    return SimpleWeightedEdge(T(t[1]), T(t[2]), U(t[3]))
end

function SimpleWeightedEdge{T,U}(t::NTuple{2}) where {T<:Integer} where {U<:Real}
    return SimpleWeightedEdge(T(t[1]), T(t[2]), one(U))
end

function SimpleWeightedEdge{T,U}(x, y) where {T<:Integer} where {U<:Real}
    return SimpleWeightedEdge(x, y, one(U))
end

SimpleWeightedEdge(x, y) = SimpleWeightedEdge(x, y, one(Float64))

Base.eltype(::AbstractSimpleWeightedEdge{T}) where {T} = T

# Accessors
Graphs.src(e::AbstractSimpleWeightedEdge) = e.src
Graphs.dst(e::AbstractSimpleWeightedEdge) = e.dst
weight(e::AbstractSimpleWeightedEdge) = e.weight

# I/O
function Base.show(io::IO, e::AbstractSimpleWeightedEdge)
    return print(io, "Edge $(e.src) => $(e.dst) with weight $(e.weight)")
end

# Conversions
Base.Tuple(e::AbstractSimpleWeightedEdge) = (src(e), dst(e), weight(e))

function (::Type{SimpleWeightedEdge{T,U}})(
    e::AbstractSimpleWeightedEdge
) where {T<:Integer,U<:Real}
    return SimpleWeightedEdge{T,U}(T(e.src), T(e.dst), U(e.weight))
end

# Convenience functions - note that these do not use weight.
Base.reverse(e::T) where {T<:AbstractSimpleWeightedEdge} = T(dst(e), src(e), weight(e))

function Base.:(==)(e1::AbstractSimpleWeightedEdge, e2::AbstractSimpleWeightedEdge)
    return (src(e1) == src(e2) && dst(e1) == dst(e2))
end

function Base.:(==)(e1::AbstractSimpleWeightedEdge, e2::AbstractEdge)
    return (src(e1) == src(e2) && dst(e1) == dst(e2))
end

function Base.:(==)(e1::AbstractEdge, e2::AbstractSimpleWeightedEdge)
    return (src(e1) == src(e2) && dst(e1) == dst(e2))
end
