@testset "A*" begin
    g = SimpleWeightedGraph(3)  # or use `SimpleWeightedDiGraph` for directed graphs
    add_edge!(g, 1, 2, 0.5)
    add_edge!(g, 2, 3, 0.8)
    add_edge!(g, 1, 3, 2.0)
    @test_broken length(a_star(g, 1, 3)) == 2
    distmx = weights(g)
    heuristic(v) = 0
    edgetype_to_return = Graphs.SimpleEdge
    @test length(a_star(g, 1, 3, distmx, heuristic, edgetype_to_return)) == 2
end
